#!/usr/bin/env python3
import os
import json
import csv
import re
import argparse
from typing import Optional, Tuple


BINDER_REGEX = re.compile(r"^AB\d+$")   # AB1, AB2, ...


def find_binder_root(start: str, max_up: int = 6) -> Optional[str]:
    """Walk up until we hit AB{number}."""
    cur = os.path.abspath(start)
    for _ in range(max_up):
        if BINDER_REGEX.match(os.path.basename(cur)):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    return None


def find_input_json(binder_root: str) -> Optional[str]:
    """Find AF3 input JSON (contains 'sequences')."""
    for fn in os.listdir(binder_root):
        if fn.endswith(".json"):
            path = os.path.join(binder_root, fn)
            try:
                with open(path) as f:
                    data = json.load(f)
                if "sequences" in data:
                    return path
            except Exception:
                pass
    return None


def extract_BC_sequences(json_path: str) -> Tuple[str, str]:
    """Extract chain B (heavy) and C (light) sequences."""
    with open(json_path) as f:
        data = json.load(f)

    heavy, light = "", ""

    for entry in data.get("sequences", []):
        prot = entry.get("protein", {})
        ids = prot.get("id", [])
        if isinstance(ids, str):
            ids = [ids]
        seq = prot.get("sequence", "")
        if "B" in ids:
            heavy = seq
        if "C" in ids:
            light = seq

    return heavy, light


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?", default=".")
    ap.add_argument("--threshold", type=float, default=0.35)
    ap.add_argument("--out", default="iptm_pass_binders.csv")
    args = ap.parse_args()

    rows = []

    for dirpath, _, filenames in os.walk(args.root):
        for fn in filenames:
            if not fn.endswith("_summary_confidences.json"):
                continue

            with open(os.path.join(dirpath, fn)) as f:
                conf = json.load(f)

            iptm = conf.get("iptm")
            if iptm is None or iptm <= args.threshold:
                continue

            binder_root = find_binder_root(dirpath)
            if binder_root is None:
                continue

            input_json = find_input_json(binder_root)
            if input_json is None:
                continue

            heavy, light = extract_BC_sequences(input_json)

            rows.append({
                "binder": os.path.basename(binder_root),
                "iptm": iptm,
                "heavy_chain_B": heavy,
                "light_chain_C": light,
            })

    if not rows:
        print("No models passed the iptm cutoff.")
        return

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["binder", "iptm", "heavy_chain_B", "light_chain_C"]
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} entries → {args.out}")


if __name__ == "__main__":
    main()

