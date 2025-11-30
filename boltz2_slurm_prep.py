#!/usr/bin/env python3
import argparse
import os
import stat


def parse_fasta(fasta_path):
    """
    Parse a FASTA file -> list of (header, sequence) tuples.
    """
    records = []
    header = None
    seq_chunks = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # save previous
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:].split()[0]  # take first token after '>'
                seq_chunks = []
            else:
                seq_chunks.append(line)

    # last record
    if header is not None:
        records.append((header, "".join(seq_chunks)))

    return records


def make_yaml_text(chain_a_seq, chain_b_seq):
    """
    Build the YAML text for Boltz-2.

    Minimal structure: just two protein chains A (fixed) and B (binder).
    We rely on --use_msa_server to generate MSAs automatically.
    """
    yaml = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {chain_a_seq}
  - protein:
      id: B
      sequence: {chain_b_seq}
"""
    return yaml


def make_sbatch_script(job_name, yaml_filename, output_dir):
    """
    Build the SBATCH script text for Boltz-2.
    Adjust module load / environment to your cluster setup if needed.
    """
    script = f"""#!/bin/bash
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

cd $SLURM_SUBMIT_DIR
mkdir -p {output_dir}

boltz predict ./{yaml_filename} \\
  --use_msa_server \\
  --out_dir ./
"""
    return script


def main():
    parser = argparse.ArgumentParser(
        description="Generate Boltz-2 YAML + SLURM scripts for a list of binders."
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Multi-FASTA file containing binder (chain B) sequences.",
    )
    parser.add_argument(
        "--chain-a-seq",
        required=True,
        help="Sequence for chain A (same for all jobs, e.g. antigen).",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Base output directory (default: current directory).",
    )
    parser.add_argument(
        "--name-prefix",
        default="binder",
        help='Prefix for job directories, e.g. "binder" -> binder_1, binder_2, ...',
    )
    args = parser.parse_args()

    binders = parse_fasta(args.fasta)
    if not binders:
        raise SystemExit("No sequences found in FASTA.")

    os.makedirs(args.outdir, exist_ok=True)

    for i, (header, seq_b) in enumerate(binders, start=1):
        dir_name = f"{args.name_prefix}_{i}"  # binder_1, binder_2, ...
        dir_path = os.path.join(args.outdir, dir_name)
        os.makedirs(dir_path, exist_ok=True)

        # YAML
        yaml_name = "input.yaml"
        yaml_path = os.path.join(dir_path, yaml_name)

        yaml_text = make_yaml_text(
            chain_a_seq=args.chain_a_seq,
            chain_b_seq=seq_b,
        )
        # Optionally add the header as a comment at the top
        yaml_text = f"# binder: {header}\n" + yaml_text

        with open(yaml_path, "w") as yf:
            yf.write(yaml_text)

        # SBATCH script
        job_name = f"BOLTZ_{dir_name}"
        sbatch_name = "submit_boltz.sbatch"
        sbatch_path = os.path.join(dir_path, sbatch_name)
        sbatch_text = make_sbatch_script(
            job_name=job_name,
            yaml_filename=yaml_name,
            output_dir="./",  # change if you want a different out_dir
        )
        with open(sbatch_path, "w") as sf:
            sf.write(sbatch_text)

        # make the submit script executable
        st = os.stat(sbatch_path)
        os.chmod(sbatch_path, st.st_mode | stat.S_IEXEC)

        print(f"Created: {dir_path} (header: {header})")

    print("Done.")


if __name__ == "__main__":
    main()
