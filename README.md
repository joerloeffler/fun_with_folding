# AF3 & Boltz Binder Design Utilities

Collection of helper scripts for setting up, running, and analyzing
**AlphaFold3** and **Boltz-2** antibody / binder design jobs on SLURM clusters - Queues and modules to be adapted!

This repository is intentionally **script-based** for now.
A unified Python package + CLI is current work in progress.

---

## Repository Scope

This repo covers three main tasks:

1. **Job preparation**
   - Generate JSON/YAML inputs
   - Generate SLURM submission scripts
2. **Post-processing**
   - Extract **ipTM** and **ipSAE**
   - Rank / filter binders
   - Recover binder sequences from input files
3. **Batch automation**
   - Designed to run over many `binder_*` or `AB*` directories

All scripts are standalone and can be mixed into existing pipelines.

---

## Directory Conventions 

The scripts assume a directory layout like:

```
project/
├── binder_1/
│   ├── input.yaml
│   ├── input.json
│   ├── submit_*.sbatch
│   └── boltz_results_input/
│       └── predictions/input/
│           ├── pae_input_model_0.npz
│           ├── input_model_0.cif
│           └── confidence_input_model_0.json
├── binder_2/
├── binder_3/
└── ...
```

or for AF3 antigen–antibody grids:

```
AB1/
AB2/
AB3/
...
```

---

## Scripts Overview

---

### AF3_slurm_prep.py

Generate **AlphaFold3 JSON inputs and SLURM scripts** for
single-chain binders against a fixed target.

**Chains**
- A → antigen / target
- B → binder

**Example**
```bash
python AF3_slurm_prep.py \
  --fasta binders.fasta \
  --chain-a-seq MPAENKKVRFENTTSDKG... \
  --outdir af3_runs \
  --name-prefix binder \
  --model-seed 1
```

Creates:
```
binder_1/
  ├── input.json
  └── submit_af3.sbatch
binder_2/
...
```

---

### AF3_AB_slurm_prep.py

Prepare **full antigen–antibody AF3 jobs** by iterating over all
heavy × light chain combinations.

**Chains**
- A → antigen
- B → heavy chain
- C → light chain

**Example**
```bash
python AF3_AB_slurm_prep.py \
  --chain-a-seq MPAENKKVRFENTTSDKG... \
  --heavy-fasta heavy.fasta \
  --light-fasta light.fasta \
  --outdir af3_ab_runs \
  --name-prefix AB
```

Produces:
```
AB1/
AB2/
AB3/
...
```

Each directory contains a valid AF3 `input.json` and `submit_af3.sbatch`.

---

### boltz2_slurm_prep.py

Prepare **Boltz-2 binder prediction jobs** with minimal YAML inputs.

**Chains**
- A → fixed target
- B → binder

**Example**
```bash
python boltz2_slurm_prep.py \
  --fasta binders.fasta \
  --chain-a-seq MPAENKKVRFENTTSDKG... \
  --outdir boltz_runs \
  --name-prefix binder
```

Creates:
```
binder_1/input.yaml
binder_1/submit_boltz.sbatch
```

---

### boltz_ipsae_boltz2.py

Collect **ipSAE / ipTM** from Boltz-2 predictions and recover
binder sequences.

- Runs `ipsae.py` automatically if needed
- Extracts max ipSAE
- Writes a ranked CSV

**Example**
```bash
python boltz_ipsae_boltz2.py \
  --root-dir boltz_runs \
  --ipsae-script ~/storage/Scripts/ipsae.py \
  --threshold 0.7
```

Output:
```
overview.csv
```

Columns:
```
binder_id, ipSAE, ipTM, sequence_B
```

---

### get_ipsae_AF3.py

Collect **ipTM scores** from AlphaFold3 runs and extract
heavy/light sequences from AF3 input JSON files.

**Example**
```bash
python get_ipsae_AF3.py \
  --threshold 0.35 \
  --out iptm_pass_binders.csv
```

Output:
```csv
binder,iptm,heavy_chain_B,light_chain_C
AB12,0.82,EVQLVESGGGLVK...,DIQMTQSPSSLS...
```

---

### AF3_AB_iptm.py

Batch-scan AF3 directories and extract **high-confidence complexes**.

- Searches recursively for `*_summary_confidences.json`
- Applies ipTM cutoff
- Recovers heavy/light sequences from AF3 inputs

**Example**
```bash
python AF3_AB_iptm.py \
  --threshold 0.4 \
  --out passing_af3_complexes.csv
```

---

## Dependencies

Minimal dependencies; most scripts rely only on the Python standard library.

Optional:
```bash
pip install pyyaml
```

External tools assumed:
- AlphaFold3
- Boltz-2
- ipsae.py

---

## Note

This repository reflects ** HPC workflows**, not a polished library.
Scripts are meant to be copied, adapted, and extended.

---

## Roadmap

Planned improvements:

- Unified CLI (`af3-boltz`)
- Remove hard coded paths (flexibility!)
- Proper Python package layout
- Visualization & ranking helpers

---


