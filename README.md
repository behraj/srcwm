# SRCWM — Shared-Rescue Cell World Model (E0 calibration)

Cell-type-resolved, multi-disease drug repurposing via action-conditioned latent
dynamics with honest zero-shot evaluation. This repository currently contains the
**E0 calibration phase** (`EXPERIMENTAL_DESIGN.md`, Milestone 0) as a single
self-contained notebook: `srcwm_e0_calibration.ipynb`.

## What the notebook does

- **Installs its own dependencies** (numpy, pandas, scipy, scikit-learn, h5py,
  matplotlib, cellxgene-census) in the first cell.
- **Downloads every dataset it uses at runtime** and freezes a snapshot:
  - LINCS L1000 Phase 1 metadata (GSE92742) — NCBI FTP, SHA512-verified against
    the published `GSE92742_SHA512SUMS.txt.gz`.
  - GEO calibration series `GSE5281` (Alzheimer's), `GSE55235` (rheumatoid
    arthritis), `GSE24206` (pulmonary fibrosis) with NCBI platform annotations.
  - STRING v12.0 human PPI (9606).
  - ChEMBL `drug_indication` + `molecule` paged REST snapshot (ground truth).
  - CZ CELLxGENE Census `2025-11-08` LTS slices (brain / lung / blood).
- **Validates the toolchain on a known-answer synthetic world** (E4 smoke):
  planted shared-rescue drugs must be recovered by the SRCWM `min_k` ranking on
  every LOGO fold; negative controls (scrambled delta, scrambled labels) must
  collapse. All gates are `assert`ed — the notebook aborts, it does not skip.
- **Does not fake anything**: deep-model arms (M-CPA/M-GEARS/M-SCGPT/M-STATE) are
  specified but deliberately not executed; failures are printed as failures.

## Two profiles

| | `SMOKE = True` (default) | `SMOKE = False` (Lambda FULL) |
|---|---|---|
| LINCS Level-5 GCTX (21 GB) | not downloaded | downloaded + SHA512-verified |
| Real-panel scoring (E4) | synthetic world only | all approved perturbagens, LOGO folds |
| Bootstrap | B = 2,000 | B = 10,000 (pre-registered E5) |
| Disk needs | < 1 GB | ~60 GB recommended |
|_runtime | ~15–25 min | ~21 GB download + several CPU hours |

Everything that executes in SMOKE is real: real GEO data, real ChEMBL API, real
Census slices, real STRING. The only synthetic components are the explicitly
labelled known-answer worlds used for toolchain calibration.

## Running

```bash
git clone https://github.com/behraj/srcwm && cd srcwm
# open srcwm_e0_calibration.ipynb in JupyterLab / Jupyter and run all cells
```

On a Lambda server: `pip install jupyterlab` then launch, or
`jupyter nbconvert --to notebook --execute --inplace srcwm_e0_calibration.ipynb`.
Set `SMOKE = False` in the config cell for the FULL profile. All outputs land in
`results/` and `figures/`, with `data/SNAPSHOT.md` (E10) recording sources,
versions, and SHA-256 checksums.

## Repository layout

```
srcwm_e0_calibration.ipynb      the notebook (run this; clean, no outputs)
srcwm_e0_calibration_evidence.ipynb  a smoke-executed copy with outputs (evidence)
cells_src/                      one file per notebook cell (buildable source)
assemble_notebook.py            rebuilds the .ipynb from cells_src/
EXPERIMENTAL_DESIGN.md          the pre-registered experimental design (E0–E13)
results figures                 produced by executing the notebook (committed
                                as execution evidence; data/ is git-ignored)
```

To regenerate the notebook from sources:
`pip install nbformat && python assemble_notebook.py`.
