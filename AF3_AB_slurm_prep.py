#!/usr/bin/env python3
import argparse
import os
import json
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


def make_json(antigen_seq, heavy_seq, light_seq, name, model_seed=1):
    """
    Build the JSON dict for AF3 for an antigen–antibody complex:
      - chain A = antigen
      - chain B = heavy chain
      - chain C = light chain
    """
    return {
        "name": name,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],  # Antigen
                    "sequence": antigen_seq,
                }
            },
            {
                "protein": {
                    "id": ["B"],  # Heavy chain
                    "sequence": heavy_seq,
                }
            },
            {
                "protein": {
                    "id": ["C"],  # Light chain
                    "sequence": light_seq,
                }
            },
        ],
        "modelSeeds": [model_seed],
        "dialect": "alphafold3",
        "version": 1,
    }


def make_sbatch_script(job_name, json_filename, output_dir):
    """
    Build the SBATCH script text.
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
module load alphafold/3.0

cd $SLURM_SUBMIT_DIR
mkdir -p {output_dir}

python `which run_alphafold.py` \\
  --json_path=./{json_filename} \\
  --output_dir=./
"""
    return script


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate AlphaFold3 JSON + SLURM scripts for all "
            "antigen–heavy–light antibody combinations."
        )
    )
    parser.add_argument(
        "--chain-a-seq",
        required=True,
        help="Antigen sequence for chain A (same for all jobs).",
    )
    parser.add_argument(
        "--heavy-fasta",
        required=True,
        help="Multi-FASTA file containing heavy chain sequences.",
    )
    parser.add_argument(
        "--light-fasta",
        required=True,
        help="Multi-FASTA file containing light chain sequences.",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Base output directory (default: current directory).",
    )
    parser.add_argument(
        "--name-prefix",
        default="AB",
        help='Prefix for job directories, e.g. "AB" -> AB1, AB2, ...',
    )
    parser.add_argument(
        "--model-seed",
        type=int,
        default=1,
        help="Model seed to use in JSON (default: 1).",
    )
    args = parser.parse_args()

    antigen_seq = args.chain_a_seq

    heavies = parse_fasta(args.heavy_fasta)
    lights = parse_fasta(args.light_fasta)

    if not heavies:
        raise SystemExit("No sequences found in HEAVY FASTA.")
    if not lights:
        raise SystemExit("No sequences found in LIGHT FASTA.")

    os.makedirs(args.outdir, exist_ok=True)

    combo_index = 1
    for h_header, h_seq in heavies:
        for l_header, l_seq in lights:
            dir_name = f"{args.name_prefix}{combo_index}"  # AB1, AB2, ...
            dir_path = os.path.join(args.outdir, dir_name)
            os.makedirs(dir_path, exist_ok=True)

            # JSON
            json_name = "input.json"
            json_path = os.path.join(dir_path, json_name)

            # Name can encode which heavy/light we used
            combo_name = f"{h_header}__{l_header}"

            json_data = make_json(
                antigen_seq=antigen_seq,
                heavy_seq=h_seq,
                light_seq=l_seq,
                name=combo_name,
                model_seed=args.model_seed,
            )
            with open(json_path, "w") as jf:
                json.dump(json_data, jf, indent=2)

            # SBATCH script
            job_name = f"AF3_{dir_name}"
            sbatch_name = "submit_af3.sbatch"
            sbatch_path = os.path.join(dir_path, sbatch_name)
            sbatch_text = make_sbatch_script(
                job_name=job_name,
                json_filename=json_name,
                output_dir="./",
            )
            with open(sbatch_path, "w") as sf:
                sf.write(sbatch_text)

            # make the submit script executable
            st = os.stat(sbatch_path)
            os.chmod(sbatch_path, st.st_mode | stat.S_IEXEC)

            print(
                f"Created: {dir_path} (Ag: chain A, heavy: {h_header} -> B, light: {l_header} -> C)"
            )

            combo_index += 1

    print("Done.")


if __name__ == "__main__":
    main()
