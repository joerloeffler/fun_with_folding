#!/usr/bin/env python3
"""
Boltz-2 pipeline generator with Slurm dependencies (MSA -> Predict) + ions-as-ligands.

Fixes included (per your requests):
1) Restraints/constraints YAML format exactly like:
   constraints:
     - contact:
         token1: [A, 218]
         token2: [B, 10]
         max_distance: 10.0
         force: true

2) --no-msa-server now changes YAML + jobdir contents:
   - Default (MSA server ON): YAML has NO "msa:" fields, NO msa.sbatch is written,
     submit_chain.sh submits predict only, predict.sbatch includes --use_msa_server.
   - With --no-msa-server: YAML includes "msa:" paths, msa.sbatch is written,
     submit_chain.sh submits msa then predict with dependency, predict.sbatch has no --use_msa_server.

Ions:
- Implemented as ligands:
    - ligand:
        id: I1
        smiles: "[Cu+2]"
"""

import argparse
import os
import stat
from typing import List, Tuple, Optional, Dict, Any


# ----------------------------
# Presets (Å)
# ----------------------------
PP_CONTACT_PRESET_MAXDIST = {"loose": 10.0, "medium": 8.0, "tight": 6.0}
POCKET_PRESET_MAXDIST = {"loose": 8.0, "medium": 6.0, "tight": 4.5}


# ----------------------------
# Parsers
# ----------------------------
def parse_fasta(fasta_path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    seq_chunks: List[str] = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        records.append((header, "".join(seq_chunks)))
    return records


def parse_smi(smi_path: str) -> List[Tuple[str, str]]:
    ligands: List[Tuple[str, str]] = []
    auto_i = 0
    with open(smi_path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0].strip()
            if not smiles:
                continue
            if len(parts) >= 2:
                name = parts[1].strip()
            else:
                auto_i += 1
                name = f"lig{auto_i}"
            ligands.append((name, smiles))
    return ligands


def parse_pp_contacts_arg(s: str) -> List[Tuple[Tuple[str, int], Tuple[str, int]]]:
    out: List[Tuple[Tuple[str, int], Tuple[str, int]]] = []
    items = [x.strip() for x in s.split(",") if x.strip()]
    for it in items:
        if "-" not in it:
            raise ValueError(f"Bad --pp-contacts token '{it}'. Expected like A:177-B:12")
        left, right = it.split("-", 1)
        if ":" not in left or ":" not in right:
            raise ValueError(f"Bad --pp-contacts token '{it}'. Expected like A:177-B:12")
        ch1, r1 = left.split(":", 1)
        ch2, r2 = right.split(":", 1)
        out.append(((ch1.strip(), int(r1)), (ch2.strip(), int(r2))))
    return out


def parse_pocket_contacts_arg(s: str) -> List[Tuple[str, int]]:
    out: List[Tuple[str, int]] = []
    items = [x.strip() for x in s.split(",") if x.strip()]
    for it in items:
        if ":" not in it:
            raise ValueError(f"Bad --pocket-contacts token '{it}'. Expected like A:177")
        ch, r = it.split(":", 1)
        out.append((ch.strip(), int(r)))
    return out


# ----------------------------
# Generic restraintlist parser (optional)
# ----------------------------
def parse_restraintlist(s: str) -> List[Dict[str, Any]]:
    """
    Semicolon-separated entries:
      contact:A:218:B:10[:maxdist][:force]
      pocket:L:A:177,A:178[:maxdist][:force]

    force accepts: true/false (case-insensitive)
    """
    out: List[Dict[str, Any]] = []
    entries = [e.strip() for e in s.split(";") if e.strip()]

    def parse_bool(x: str) -> bool:
        xl = x.strip().lower()
        if xl in ("true", "t", "1", "yes", "y"):
            return True
        if xl in ("false", "f", "0", "no", "n"):
            return False
        raise ValueError(f"Bad boolean '{x}' (use true/false).")

    for e in entries:
        parts = [p.strip() for p in e.split(":") if p.strip()]
        kind = parts[0].lower()

        if kind == "contact":
            if len(parts) < 5:
                raise ValueError(
                    f"Bad restraint entry '{e}'. Expected contact:A:218:B:10[:maxdist][:force]"
                )
            ch1, r1, ch2, r2 = parts[1], int(parts[2]), parts[3], int(parts[4])
            maxdist = float(parts[5]) if len(parts) >= 6 else None
            force = parse_bool(parts[6]) if len(parts) >= 7 else None
            out.append(
                {
                    "contact": {
                        "token1": [ch1, r1],
                        "token2": [ch2, r2],
                        "max_distance": maxdist,
                        "force": force,
                    }
                }
            )

        elif kind == "pocket":
            if len(parts) < 3:
                raise ValueError(
                    f"Bad restraint entry '{e}'. Expected pocket:L:A:177,A:178[:maxdist][:force]"
                )
            binder = parts[1]
            contacts_str = parts[2]
            contacts: List[List[Any]] = []
            for tok in [t.strip() for t in contacts_str.split(",") if t.strip()]:
                if ":" not in tok:
                    raise ValueError(f"Bad pocket contact '{tok}' in '{e}' (expected A:177).")
                ch, rr = tok.split(":", 1)
                contacts.append([ch.strip(), int(rr)])
            maxdist = float(parts[3]) if len(parts) >= 4 else None
            force = parse_bool(parts[4]) if len(parts) >= 5 else None
            out.append(
                {
                    "pocket": {
                        "binder": binder,
                        "contacts": contacts,
                        "max_distance": maxdist,
                        "force": force,
                    }
                }
            )

        else:
            raise ValueError(f"Unknown restraint kind '{parts[0]}' in '{e}'. Supported: contact, pocket.")

    return out


# ----------------------------
# YAML helpers/builders
# ----------------------------
def yaml_escape_single_quotes(s: str) -> str:
    return s.replace("'", "''")


def yaml_bool(b: bool) -> str:
    return "true" if b else "false"


def ion_smiles(ion_type: str, ion_charge: int) -> str:
    """
    Build bracketed ion SMILES like:
      CU +2 -> [Cu+2]
      CL -1 -> [Cl-]
      NA +1 -> [Na+]
    """
    sym = ion_type.strip()
    if not sym:
        raise ValueError("ion_type is empty")
    sym_norm = sym[0].upper() + sym[1:].lower() if len(sym) > 1 else sym.upper()

    q = int(ion_charge)
    if q == 0:
        raise ValueError("ion_charge cannot be 0 for an ion")

    if abs(q) == 1:
        sign = "+" if q > 0 else "-"
        charge_txt = sign
    else:
        sign = "+" if q > 0 else "-"
        charge_txt = f"{sign}{abs(q)}"

    return f"[{sym_norm}{charge_txt}]"


def ion_ligands_to_yaml_lines(ion_type: str, ion_charge: int, ion_quant: int, start_index: int = 1) -> List[str]:
    """
    Emit N ligand blocks for ions:
      - ligand:
          id: I1
          smiles: "[Cu+2]"
    """
    n = int(ion_quant)
    if n <= 0:
        raise ValueError("ion_quant must be >= 1")

    smi = ion_smiles(ion_type, ion_charge)
    safe = yaml_escape_single_quotes(smi)

    lines: List[str] = []
    for k in range(start_index, start_index + n):
        lines.append("  - ligand:")
        lines.append(f"      id: I{k}")
        lines.append(f"      smiles: '{safe}'")
    return lines


def constraints_to_yaml_lines(
    constraints: List[Dict[str, Any]],
    *,
    default_contact_maxdist: float,
    default_contact_force: bool,
    default_pocket_maxdist: float,
    default_pocket_force: bool,
) -> List[str]:
    """
    Writes constraints in the exact indentation/format you asked for.
    """
    lines: List[str] = []
    lines.append("constraints:")
    for c in constraints:
        if "contact" in c:
            d = c["contact"]
            maxd = d.get("max_distance")
            if maxd is None:
                maxd = default_contact_maxdist
            force = d.get("force")
            if force is None:
                force = default_contact_force

            lines.append("  - contact:")
            lines.append(f"      token1: [{d['token1'][0]}, {int(d['token1'][1])}]")
            lines.append(f"      token2: [{d['token2'][0]}, {int(d['token2'][1])}]")
            lines.append(f"      max_distance: {float(maxd)}")
            lines.append(f"      force: {yaml_bool(bool(force))}")

        elif "pocket" in c:
            d = c["pocket"]
            maxd = d.get("max_distance")
            if maxd is None:
                maxd = default_pocket_maxdist
            force = d.get("force")
            if force is None:
                force = default_pocket_force

            lines.append("  - pocket:")
            lines.append(f"      binder: {d['binder']}")
            lines.append("      contacts:")
            for ch, rr in d["contacts"]:
                lines.append(f"        - [{ch}, {int(rr)}]")
            lines.append(f"      max_distance: {float(maxd)}")
            lines.append(f"      force: {yaml_bool(bool(force))}")

        else:
            raise ValueError(f"Unknown constraint dict: {c}")

    return lines


def make_yaml_text_binder(
    chain_a_seq: str,
    chain_b_seq: str,
    binder_header: str,
    msa_a: Optional[str],
    msa_b: Optional[str],
    pp_contacts: Optional[List[Tuple[Tuple[str, int], Tuple[str, int]]]] = None,
    pp_max_distance: float = 6.0,
    pp_force: bool = False,
    extra_constraints: Optional[List[Dict[str, Any]]] = None,
    ions: bool = False,
    ion_type: Optional[str] = None,
    ion_charge: Optional[int] = None,
    ion_quant: Optional[int] = None,
) -> str:
    lines: List[str] = []
    lines.append(f"# binder: {binder_header}")
    lines.append("version: 1")
    lines.append("sequences:")

    # proteins
    lines.append("  - protein:")
    lines.append("      id: A")
    if msa_a:
        lines.append(f"      msa: {msa_a}")
    lines.append(f"      sequence: {chain_a_seq}")

    lines.append("  - protein:")
    lines.append("      id: B")
    if msa_b:
        lines.append(f"      msa: {msa_b}")
    lines.append(f"      sequence: {chain_b_seq}")

    # ions as ligands
    if ions:
        if ion_type is None or ion_charge is None or ion_quant is None:
            raise ValueError("ions=True requires ion_type, ion_charge, ion_quant")
        lines.extend(ion_ligands_to_yaml_lines(ion_type, ion_charge, ion_quant, start_index=1))

    # constraints
    constraints: List[Dict[str, Any]] = []

    if pp_contacts:
        for (t1, t2) in pp_contacts:
            ch1, r1 = t1
            ch2, r2 = t2
            constraints.append(
                {
                    "contact": {
                        "token1": [ch1, int(r1)],
                        "token2": [ch2, int(r2)],
                        "max_distance": float(pp_max_distance),
                        "force": bool(pp_force),
                    }
                }
            )

    if extra_constraints:
        constraints.extend(extra_constraints)

    if constraints:
        lines.extend(
            constraints_to_yaml_lines(
                constraints,
                default_contact_maxdist=pp_max_distance,
                default_contact_force=pp_force,
                default_pocket_maxdist=POCKET_PRESET_MAXDIST["tight"],
                default_pocket_force=False,
            )
        )
    else:
        lines.append("constraints: []")

    return "\n".join(lines) + "\n"


def make_yaml_text_ligand(
    chain_a_seq: str,
    ligand_name: str,
    ligand_smiles: str,
    msa_a: Optional[str],
    pocket_contacts: Optional[List[Tuple[str, int]]] = None,
    pocket_max_distance: float = 4.5,
    pocket_force: bool = False,
    extra_constraints: Optional[List[Dict[str, Any]]] = None,
    ions: bool = False,
    ion_type: Optional[str] = None,
    ion_charge: Optional[int] = None,
    ion_quant: Optional[int] = None,
) -> str:
    safe_smiles = yaml_escape_single_quotes(ligand_smiles)

    lines: List[str] = []
    lines.append(f"# ligand: {ligand_name} | smiles: {ligand_smiles}")
    lines.append("version: 1")
    lines.append("sequences:")

    # protein
    lines.append("  - protein:")
    lines.append("      id: A")
    if msa_a:
        lines.append(f"      msa: {msa_a}")
    lines.append(f"      sequence: {chain_a_seq}")

    # main ligand
    lines.append("  - ligand:")
    lines.append("      id: L")
    lines.append(f"      smiles: '{safe_smiles}'")

    # ions as ligands
    if ions:
        if ion_type is None or ion_charge is None or ion_quant is None:
            raise ValueError("ions=True requires ion_type, ion_charge, ion_quant")
        lines.extend(ion_ligands_to_yaml_lines(ion_type, ion_charge, ion_quant, start_index=1))

    # constraints
    constraints: List[Dict[str, Any]] = []

    if pocket_contacts:
        contacts = [[ch, int(r)] for ch, r in pocket_contacts]
        constraints.append(
            {
                "pocket": {
                    "binder": "L",
                    "contacts": contacts,
                    "max_distance": float(pocket_max_distance),
                    "force": bool(pocket_force),
                }
            }
        )

    if extra_constraints:
        constraints.extend(extra_constraints)

    if constraints:
        lines.extend(
            constraints_to_yaml_lines(
                constraints,
                default_contact_maxdist=PP_CONTACT_PRESET_MAXDIST["tight"],
                default_contact_force=False,
                default_pocket_maxdist=pocket_max_distance,
                default_pocket_force=pocket_force,
            )
        )
    else:
        lines.append("constraints: []")

    return "\n".join(lines) + "\n"


# ----------------------------
# File builders (MSA sbatch + Predict sbatch + submit scripts)
# ----------------------------
def make_protein_chains_fasta(chain_a_seq: str, chain_b_seq: Optional[str]) -> str:
    lines = [">A", chain_a_seq]
    if chain_b_seq is not None:
        lines += [">B", chain_b_seq]
    return "\n".join(lines) + "\n"


def make_msa_sbatch_script(
    *,
    job_name: str,
    micromamba_hook_cmd: str,
    micromamba_env: str,
    boltz_nvs_py: str,
    protein_chains: str,
    partition: str,
    nodelist: str,
    nodes: int,
    cpus_per_task: int,
    mem: str,
) -> str:
    nodelist_line = f"#SBATCH --nodelist={nodelist}\n" if nodelist.strip() else ""
    return f"""#!/bin/bash
#SBATCH --partition=ampere
#SBATCH --ntasks-per-node=8
#SBATCH --mem=25gb
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --job-name={job_name}
#SBATCH --exclude=nodea0401,nodea0403

module purge
# Adjust this to whatever module/conda env you use for Boltz-2
module load boltz

python {boltz_nvs_py} --faa ./protein_chains.fasta --outdir ./boltz2output --msadir ./boltz2output/msas --protein-chains={protein_chains}
"""


def make_predict_sbatch_script(
    *,
    job_name: str,
    micromamba_hook_cmd: str,
    micromamba_env: str,
    boltz_bin: str,
    yaml_filename: str,
    out_dir: str,
    use_msa_server: bool,
    partition: str,
    nodelist: str,
    nodes: int,
    cpus_per_task: int,
    mem: str,
    gres: str,
) -> str:
    args = [boltz_bin, "predict", f"./{yaml_filename}"]
    if use_msa_server:
        args.append("--use_msa_server")
    args += ["--out_dir", out_dir]
    pretty = " \\\n  ".join(args)

    nodelist_line = f"#SBATCH --nodelist={nodelist}\n" if nodelist.strip() else ""

    return f"""#!/bin/bash
#SBATCH --partition=ampere
#SBATCH --ntasks-per-node=8
#SBATCH --mem=25gb
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --job-name={job_name}
#SBATCH --exclude=nodea0401,nodea0403

module purge
# Adjust this to whatever module/conda env you use for Boltz-2
module load boltz

{pretty}
"""


def make_submit_chain_sh(msa_sbatch: str = "msa.sbatch", pred_sbatch: str = "predict.sbatch") -> str:
    return f"""#!/bin/bash
set -euo pipefail

jid=$(sbatch --parsable {msa_sbatch})
sbatch --dependency=afterok:$jid {pred_sbatch}

echo "Submitted dependency chain:"
echo "  MSA jobid:    $jid"
echo "  Predict job:  depends on afterok:$jid"
"""


def make_submit_predict_only_sh(pred_sbatch: str = "predict.sbatch") -> str:
    return f"""#!/bin/bash
set -euo pipefail

sbatch {pred_sbatch}
echo "Submitted predict only (MSA server mode)."
"""


def make_submit_all_sh(outdir_abs: str) -> str:
    return f"""#!/bin/bash
set -euo pipefail

base="{outdir_abs}"

for d in "$base"/*; do
  if [[ -d "$d" && -x "$d/submit_chain.sh" ]]; then
    echo "[SUBMIT] $d"
    ( cd "$d" && ./submit_chain.sh )
  fi
done

echo "Done submitting all job dirs."
"""


# ----------------------------
# IO helpers
# ----------------------------
def write_text(path: str, text: str, make_exec: bool = False) -> None:
    with open(path, "w") as f:
        f.write(text)
    if make_exec:
        st = os.stat(path)
        os.chmod(path, st.st_mode | stat.S_IEXEC)


# ----------------------------
# Main
# ----------------------------
def main():
    p = argparse.ArgumentParser(
        description="Generate Boltz-2 job dirs with msa.sbatch -> predict.sbatch dependencies, ions-as-ligands, per-dir submit_chain.sh and optional top-level submit_all.sh."
    )

    p.add_argument("--chain-a-seq", required=True)
    p.add_argument("--outdir", default=".")
    p.add_argument("--name-prefix", default="job")

    # Mode
    p.add_argument("--fasta", default=None)
    p.add_argument("--smi", default=None)

    # Distance presets + simple restraints
    p.add_argument("--preset", choices=["loose", "medium", "tight"], default="tight")
    p.add_argument("--pp-contacts", default=None)
    p.add_argument("--pp-force", action="store_true")
    p.add_argument("--pp-max-distance", type=float, default=None)

    p.add_argument("--pocket-contacts", default=None)
    p.add_argument("--pocket-force", action="store_true")
    p.add_argument("--pocket-max-distance", type=float, default=None)

    # Generic extra restraints
    p.add_argument("--restraints", action="store_true")
    p.add_argument("--restraintlist", default="")

    # Ions-as-ligands
    p.add_argument("--ions", action="store_true")
    p.add_argument("--ion_type", default=None)
    p.add_argument("--ion_charge", type=int, default=None)
    p.add_argument("--ion_quant", type=int, default=None)

    # Toolchain
    p.add_argument("--micromamba-hook", default="/z/linux/bin/micromamba shell hook -s posix")
    p.add_argument("--micromamba-env", default="/z/bio/biotools/miniconda/envs/boltz2/")
    p.add_argument("--boltz-nvs", default="/z/bio/biotools/boltz/boltz_nvs.py")
    p.add_argument("--boltz-bin", default="boltz")

    # MSA SBATCH resources (defaults: partition high)
    p.add_argument("--msa-jobname", default="MSA")
    p.add_argument("--msa-partition", default="high")
    p.add_argument("--msa-nodelist", default="")  # leave empty by default
    p.add_argument("--msa-nodes", type=int, default=1)
    p.add_argument("--msa-cpus", type=int, default=16)
    p.add_argument("--msa-mem", default="25gb")

    # Predict SBATCH resources (defaults: partition high)
    p.add_argument("--pred-jobname", default="Boltz2")
    p.add_argument("--pred-partition", default="high")
    p.add_argument("--pred-nodelist", default="ai,pika")
    p.add_argument("--pred-nodes", type=int, default=1)
    p.add_argument("--pred-cpus", type=int, default=16)
    p.add_argument("--pred-mem", default="1gb")
    p.add_argument("--pred-gres", default="gpu:1")

    # Predict options
    p.add_argument("--no-msa-server", action="store_true")
    p.add_argument("--boltz-outdir", default="./")

    # Submission helpers
    p.add_argument("--write-submit-all", action="store_true")

    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    pp_maxd = args.pp_max_distance if args.pp_max_distance is not None else PP_CONTACT_PRESET_MAXDIST[args.preset]
    pocket_maxd = args.pocket_max_distance if args.pocket_max_distance is not None else POCKET_PRESET_MAXDIST[args.preset]

    extra_constraints: Optional[List[Dict[str, Any]]] = None
    if args.restraints:
        if not args.restraintlist.strip():
            raise SystemExit("--restraints requires --restraintlist")
        extra_constraints = parse_restraintlist(args.restraintlist)

    if args.ions and (args.ion_type is None or args.ion_charge is None or args.ion_quant is None):
        raise SystemExit("--ions requires --ion_type, --ion_charge, --ion_quant")

    # If MSA server is used: omit msa paths in YAML and skip msa.sbatch entirely
    use_msa_server = not args.no_msa_server
    write_msa_job = not use_msa_server

    if use_msa_server:
        msa_a = None
        msa_b = None
    else:
        msa_a = "./boltz2output/msas/A.a3m"
        msa_b = "./boltz2output/msas/B.a3m"

    def write_job_dir(
        dir_path: str,
        yaml_text: str,
        chain_b_seq: Optional[str],
        protein_chains: str,
        *,
        write_msa_job: bool,
    ):
        os.makedirs(dir_path, exist_ok=True)

        write_text(
            os.path.join(dir_path, "protein_chains.fasta"),
            make_protein_chains_fasta(args.chain_a_seq, chain_b_seq),
        )

        # Only write msa.sbatch when NOT using the MSA server
        if write_msa_job:
            write_text(
                os.path.join(dir_path, "msa.sbatch"),
                make_msa_sbatch_script(
                    job_name=args.msa_jobname,
                    micromamba_hook_cmd=args.micromamba_hook,
                    micromamba_env=args.micromamba_env,
                    boltz_nvs_py=args.boltz_nvs,
                    protein_chains=protein_chains,
                    partition=args.msa_partition,
                    nodelist=args.msa_nodelist,
                    nodes=args.msa_nodes,
                    cpus_per_task=args.msa_cpus,
                    mem=args.msa_mem,
                ),
                make_exec=True,
            )

        write_text(os.path.join(dir_path, "input.yaml"), yaml_text)

        write_text(
            os.path.join(dir_path, "predict.sbatch"),
            make_predict_sbatch_script(
                job_name=args.pred_jobname,
                micromamba_hook_cmd=args.micromamba_hook,
                micromamba_env=args.micromamba_env,
                boltz_bin=args.boltz_bin,
                yaml_filename="input.yaml",
                out_dir=args.boltz_outdir,
                use_msa_server=use_msa_server,
                partition=args.pred_partition,
                nodelist=args.pred_nodelist,
                nodes=args.pred_nodes,
                cpus_per_task=args.pred_cpus,
                mem=args.pred_mem,
                gres=args.pred_gres,
            ),
            make_exec=True,
        )

        # submit_chain.sh depends on mode
        if write_msa_job:
            submit_script = make_submit_chain_sh()
        else:
            submit_script = make_submit_predict_only_sh()

        write_text(os.path.join(dir_path, "submit_chain.sh"), submit_script, make_exec=True)

        print(f"Created: {dir_path}")

    # Ligand mode
    if args.smi:
        ligands = parse_smi(args.smi)
        if not ligands:
            raise SystemExit(f"No ligands found in {args.smi}")

        pocket_contacts = parse_pocket_contacts_arg(args.pocket_contacts) if args.pocket_contacts else None

        for j, (lname, smiles) in enumerate(ligands, start=1):
            dir_path = os.path.join(args.outdir, f"{args.name_prefix}_lig_{j}")
            yaml_text = make_yaml_text_ligand(
                chain_a_seq=args.chain_a_seq,
                ligand_name=lname,
                ligand_smiles=smiles,
                msa_a=msa_a,
                pocket_contacts=pocket_contacts,
                pocket_max_distance=pocket_maxd,
                pocket_force=args.pocket_force,
                extra_constraints=extra_constraints,
                ions=args.ions,
                ion_type=args.ion_type,
                ion_charge=args.ion_charge,
                ion_quant=args.ion_quant,
            )
            write_job_dir(dir_path, yaml_text, chain_b_seq=None, protein_chains="A", write_msa_job=write_msa_job)

        print("Done (ligand mode).")

    # Binder mode
    else:
        if not args.fasta:
            raise SystemExit("Provide --fasta (binder mode) or --smi (ligand mode).")

        binders = parse_fasta(args.fasta)
        if not binders:
            raise SystemExit(f"No sequences found in {args.fasta}")

        pp_contacts = parse_pp_contacts_arg(args.pp_contacts) if args.pp_contacts else None

        for i, (header, seq_b) in enumerate(binders, start=1):
            dir_path = os.path.join(args.outdir, f"{args.name_prefix}_binder_{i}")
            yaml_text = make_yaml_text_binder(
                chain_a_seq=args.chain_a_seq,
                chain_b_seq=seq_b,
                binder_header=header,
                msa_a=msa_a,
                msa_b=msa_b,
                pp_contacts=pp_contacts,
                pp_max_distance=pp_maxd,
                pp_force=args.pp_force,
                extra_constraints=extra_constraints,
                ions=args.ions,
                ion_type=args.ion_type,
                ion_charge=args.ion_charge,
                ion_quant=args.ion_quant,
            )
            write_job_dir(dir_path, yaml_text, chain_b_seq=seq_b, protein_chains="A,B", write_msa_job=write_msa_job)

        print("Done (binder mode).")

    if args.write_submit_all:
        submit_all_path = os.path.join(args.outdir, "submit_all.sh")
        write_text(submit_all_path, make_submit_all_sh(os.path.abspath(args.outdir)), make_exec=True)
        print(f"Wrote: {submit_all_path}")


if __name__ == "__main__":
    main()

