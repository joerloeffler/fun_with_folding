#!/usr/bin/env python3
import os
import subprocess
import sys
import json
from typing import List, Tuple, Optional


def natural_sort_key(s: str):
    """Sort binder_1, binder_2, ..., binder_10 in human order."""
    import re
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r"(\d+)", s)]


def run_ipsae(working_dir: str, ipsae_script: str) -> None:
    """
    Run ipsae.py in the given working directory if output file is missing.
    """
    txt_name = "t0.3_seed-1_sample-0_model_10_10.txt"
    txt_path = os.path.join(working_dir, txt_name)

    if os.path.exists(txt_path):
        # Already done, skip recomputation
        print(f"[SKIP] {working_dir}: {txt_name} already exists.")
        return

    json_name = "t0.3_seed-1_sample-0_confidences.json"
    cif_name = "t0.3_seed-1_sample-0_model.cif"

    json_path = os.path.join(working_dir, json_name)
    cif_path = os.path.join(working_dir, cif_name)

    if not os.path.exists(json_path) or not os.path.exists(cif_path):
        print(f"[WARN] Missing input files in {working_dir}, skipping.")
        return

    cmd = [
        sys.executable,  # uses the same python that runs this script
        ipsae_script,
        json_name,
        cif_name,
        "10",
        "10",
    ]
    print(f"[RUN ] {' '.join(cmd)}  (cwd={working_dir})")

    try:
        subprocess.run(cmd, cwd=working_dir, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] ipsae.py failed in {working_dir}: {e}")


def parse_ipsae_and_iptm(txt_path: str) -> Tuple[float, float]:
    """
    Parse the ipsae output file and return (ipSAE, ipTM_af).

    Priority:
    1. If there is a row with Type == 'max', use that row.
    2. Otherwise, use the row with maximum ipSAE.
    """
    with open(txt_path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if len(lines) < 2:
        raise ValueError(f"No data lines in {txt_path}")

    header = lines[0].split()
    try:
        ipsae_idx = header.index("ipSAE")
        iptm_idx = header.index("ipTM_af")
        type_idx = header.index("Type")
    except ValueError as e:
        raise ValueError(f"Required column missing in {txt_path}: {e}")

    # --- First: try to find the 'max' row ---
    for line in lines[1:]:
        cols = line.split()
        if len(cols) <= max(ipsae_idx, iptm_idx, type_idx):
            continue
        if cols[type_idx] == "max":
            try:
                ipsae = float(cols[ipsae_idx])
                iptm = float(cols[iptm_idx])
            except ValueError:
                continue
            return ipsae, iptm

    # --- Fallback: use overall max ipSAE ---
    best_ipsae = None
    best_iptm = None
    for line in lines[1:]:
        cols = line.split()
        if len(cols) <= max(ipsae_idx, iptm_idx):
            continue
        try:
            val_ipsae = float(cols[ipsae_idx])
            val_iptm = float(cols[iptm_idx])
        except ValueError:
            continue
        if (best_ipsae is None) or (val_ipsae > best_ipsae):
            best_ipsae = val_ipsae
            best_iptm = val_iptm

    if best_ipsae is None or best_iptm is None:
        raise ValueError(f"No valid ipSAE/ipTM_af values in {txt_path}")

    return best_ipsae, best_iptm


def find_chain_b_sequence(binder_path: str) -> Optional[str]:
    """
    Look for a JSON file in the binder's main directory, parse it,
    and return the sequence for chain B if found.
    """
    json_candidates = [
        fname for fname in os.listdir(binder_path)
        if fname.endswith(".json") and os.path.isfile(os.path.join(binder_path, fname))
    ]

    for fname in json_candidates:
        fpath = os.path.join(binder_path, fname)
        try:
            with open(fpath, "r") as f:
                data = json.load(f)
        except Exception:
            continue

        if not isinstance(data, dict):
            continue

        sequences = data.get("sequences")
        if not isinstance(sequences, list):
            continue

        for entry in sequences:
            protein = entry.get("protein", {})
            ids = protein.get("id", [])
            if isinstance(ids, list) and "B" in ids:
                seq = protein.get("sequence")
                if isinstance(seq, str):
                    return seq

    return None


def main():
    base_dir = os.getcwd()
    ipsae_script = os.path.expanduser("~/storage/Scripts/ipsae.py")
    if not os.path.exists(ipsae_script):
        print(f"[FATAL] Cannot find ipsae.py at {ipsae_script}")
        sys.exit(1)

    binder_dirs = [
        d for d in os.listdir(base_dir)
        if d.startswith("binder_") and os.path.isdir(d)
    ]
    binder_dirs.sort(key=natural_sort_key)

    # (binder_name, ipSAE, ipTM_af, chainB_sequence or None)
    results: List[Tuple[str, float, float, Optional[str]]] = []

    for binder in binder_dirs:
        binder_path = os.path.join(base_dir, binder)
        work_dir = os.path.join(binder_path, "t0.3", "seed-1_sample-0")

        if not os.path.isdir(work_dir):
            print(f"[SKIP] {binder}: no t0.3/seed-1_sample-0 directory (not finished?).")
            continue

        # 1) Run ipsae if needed
        run_ipsae(work_dir, ipsae_script)

        # 2) Parse ipSAE & ipTM_af
        txt_path = os.path.join(work_dir, "t0.3_seed-1_sample-0_model_10_10.txt")
        if not os.path.exists(txt_path):
            print(f"[WARN] {binder}: output file {txt_path} not found, skipping.")
            continue

        try:
            ipsae, iptm = parse_ipsae_and_iptm(txt_path)
            print(f"[OK  ] {binder}: ipSAE = {ipsae:.6f}, ipTM_af = {iptm:.6f}")
        except Exception as e:
            print(f"[WARN] {binder}: failed to parse {txt_path}: {e}")
            continue

        # 3) If ipSAE is high, try to grab chain B sequence
        seq_b = None
        if ipsae > 0.75:
            seq_b = find_chain_b_sequence(binder_path)
            if seq_b is None:
                print(f"[NOTE] {binder}: ipSAE>0.75 but chain B sequence not found.")

        results.append((binder, ipsae, iptm, seq_b))

    # 4) Write overview TSV (all successful binders)
    overview_path = os.path.join(base_dir, "overview_ipsae.tsv")
    with open(overview_path, "w") as out_f:
        out_f.write("binder\tipSAE\tipTM_af\n")
        for binder, ipsae, iptm, _ in results:
            out_f.write(f"{binder}\t{ipsae:.6f}\t{iptm:.6f}\n")

    # 5) Write high-ipSAE sequences TSV
    high_path = os.path.join(base_dir, "high_ipsae_sequences.tsv")
    with open(high_path, "w") as out_f:
        out_f.write("binder\tipSAE\tipTM_af\tchainB_sequence\n")
        for binder, ipsae, iptm, seq_b in results:
            if ipsae > 0.75:
                seq_str = seq_b if seq_b is not None else ""
                out_f.write(f"{binder}\t{ipsae:.6f}\t{iptm:.6f}\t{seq_str}\n")

    print(f"\n[DONE] Wrote {len(results)} entries to {overview_path}")
    print(f"[DONE] Wrote high-ipSAE binders to {high_path}")


if __name__ == "__main__":
    main()
