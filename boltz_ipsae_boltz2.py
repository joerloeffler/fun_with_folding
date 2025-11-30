#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Optional, Tuple

try:
    import yaml  # PyYAML
except ImportError:
    print("This script requires PyYAML. Install it with: pip install pyyaml")
    sys.exit(1)


def run_ipsae(
    pred_dir: Path,
    ipsae_script: Path,
    pae_name: str = "pae_input_model_0.npz",
    cif_name: str = "input_model_0.cif",
    d0_chn: int = 10,
    d0_dom: int = 10,
) -> Optional[Path]:
    """
    Run ipsae.py in pred_dir if needed and return path to the *_d0chn_d0dom.txt output file.
    """
    pae_path = pred_dir / pae_name
    cif_path = pred_dir / cif_name

    if not pae_path.is_file() or not cif_path.is_file():
        print(f"[WARN] Missing PAE or CIF in {pred_dir}, skipping.")
        return None

    base = cif_path.stem  # e.g. "input_model_0"
    out_txt = pred_dir / f"{base}_{d0_chn}_{d0_dom}.txt"

    if not out_txt.is_file():
        print(f"[INFO] Running ipsae.py in {pred_dir}...")
        cmd = [
            "python",
            str(ipsae_script),
            str(pae_path),
            str(cif_path),
            str(d0_chn),
            str(d0_dom),
        ]
        try:
            subprocess.run(cmd, cwd=pred_dir, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] ipsae.py failed in {pred_dir}: {e}")
            return None

        if not out_txt.is_file():
            print(f"[ERROR] Expected output {out_txt} not found after ipsae.py.")
            return None

    return out_txt


def parse_ipsae_max(out_txt: Path) -> Optional[float]:
    """
    Parse ipSAE from the 'max' line in the ipsae output text file.
    """
    with out_txt.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Chn1"):
                continue
            cols = line.split()
            if len(cols) < 6:
                continue
            # Type column is index 4; ipSAE is index 5
            if cols[4] == "max":
                try:
                    return float(cols[5])
                except ValueError:
                    pass
    print(f"[WARN] No 'max' line with valid ipSAE found in {out_txt}")
    return None


def parse_iptm(pred_dir: Path, conf_name: str = "confidence_input_model_0.json") -> Optional[float]:
    """
    Parse iptm from the confidence JSON file.
    """
    conf_path = pred_dir / conf_name
    if not conf_path.is_file():
        print(f"[WARN] Missing confidence JSON in {pred_dir}")
        return None
    try:
        with conf_path.open("r") as f:
            data = json.load(f)
        # Prefer 'iptm', fall back to 'protein_iptm' if necessary
        if "iptm" in data:
            return float(data["iptm"])
        elif "protein_iptm" in data:
            return float(data["protein_iptm"])
        else:
            print(f"[WARN] No 'iptm' or 'protein_iptm' in {conf_path}")
            return None
    except Exception as e:
        print(f"[WARN] Could not parse {conf_path}: {e}")
        return None


def parse_chain_b_sequence(binder_dir: Path, yaml_name: str = "input.yaml") -> Optional[str]:
    """
    Parse chain B sequence from binder_dir/input.yaml.

    Expects structure like:
    sequences:
      - protein:
          id: A
          sequence: ...
      - protein:
          id: B
          sequence: ...
    """
    yaml_path = binder_dir / yaml_name
    if not yaml_path.is_file():
        print(f"[WARN] input.yaml not found in {binder_dir}")
        return None

    try:
        with yaml_path.open("r") as f:
            data = yaml.safe_load(f)
    except Exception as e:
        print(f"[WARN] Could not parse YAML {yaml_path}: {e}")
        return None

    sequences = data.get("sequences", [])
    for entry in sequences:
        protein = entry.get("protein", {})
        pid = protein.get("id")

        # id might be a string "B" or a list ["B"]
        if isinstance(pid, list):
            ids = [str(x) for x in pid]
        else:
            ids = [str(pid)] if pid is not None else []

        if "B" in ids:
            seq = protein.get("sequence")
            if isinstance(seq, str):
                return seq

    print(f"[WARN] No chain B sequence found in {yaml_path}")
    return None


def natural_binder_sort_key(name: str):
    """
    Sort binder directories by trailing number, e.g. binder_1, binder_2, ...
    """
    m = re.search(r"(\d+)$", name)
    if m:
        return int(m.group(1))
    return name


def collect_ipsae(
    root_dir: Path,
    ipsae_script: Path,
    threshold: float = 0.7,
    out_csv: Path = None,
) -> None:
    if out_csv is None:
        out_csv = root_dir / "overview.csv"

    rows = []

    # Find binder_* directories
    binder_dirs = [
        d for d in root_dir.iterdir()
        if d.is_dir() and d.name.startswith("binder_")
    ]
    binder_dirs = sorted(binder_dirs, key=lambda d: natural_binder_sort_key(d.name))

    print(f"[INFO] Found {len(binder_dirs)} binder directories under {root_dir}")

    for binder_dir in binder_dirs:
        binder_id = binder_dir.name
        pred_dir = binder_dir / "boltz_results_input" / "predictions" / "input"

        if not pred_dir.is_dir():
            print(f"[WARN] Prediction directory missing for {binder_id}: {pred_dir}")
            continue

        out_txt = run_ipsae(pred_dir, ipsae_script)
        if out_txt is None:
            continue

        ipsae_val = parse_ipsae_max(out_txt)
        if ipsae_val is None:
            continue

        iptm_val = parse_iptm(pred_dir)
        seq_b = parse_chain_b_sequence(binder_dir)

        print(
            f"[INFO] {binder_id}: ipSAE={ipsae_val:.4f}, "
            f"ipTM={iptm_val if iptm_val is not None else 'NA'}, "
            f"seq_B_len={len(seq_b) if seq_b else 'NA'}"
        )

        if ipsae_val >= threshold:
            rows.append({
                "binder_id": binder_id,
                "ipSAE": ipsae_val,
                "ipTM": iptm_val,
                "sequence_B": seq_b,
            })

    # Write overview.csv
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["binder_id", "ipSAE", "ipTM", "sequence_B"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"[INFO] Wrote {len(rows)} entries with ipSAE >= {threshold} to {out_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Collect ipSAE, ipTM, and chain B sequence for Boltz-2 binders."
    )
    parser.add_argument(
        "--root-dir",
        type=str,
        default=".",
        help="Root directory containing binder_* subdirectories (default: current directory).",
    )
    parser.add_argument(
        "--ipsae-script",
        type=str,
        default="~/storage/Scripts/ipsae.py",
        help="Path to ipsae.py script (default: ~/storage/Scripts/ipsae.py).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.7,
        help="ipSAE threshold for inclusion in overview.csv (default: 0.7).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional path to output CSV (default: <root-dir>/overview.csv).",
    )

    args = parser.parse_args()

    root_dir = Path(args.root_dir).resolve()
    ipsae_script = Path(os.path.expanduser(args.ipsae_script)).resolve()
    out_csv = Path(args.output).resolve() if args.output else None

    if not ipsae_script.is_file():
        print(f"[ERROR] ipsae.py not found at {ipsae_script}")
        raise SystemExit(1)

    collect_ipsae(root_dir, ipsae_script, args.threshold, out_csv)


if __name__ == "__main__":
    main()
