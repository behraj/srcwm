# SRCWM — E0 Calibration (Milestone 0)

Shared-Rescue Cell World Model (SRCWM): cell-type-resolved, multi-disease drug
repurposing via action-conditioned latent dynamics with honest zero-shot evaluation.

This notebook implements the **E0 calibration phase** of the SRCWM experimental design
(`EXPERIMENTAL_DESIGN.md` in this repository). It is a **toolchain validation**: it must
run end-to-end on a fresh machine (e.g. a Lambda GPU server), download every dataset it
consumes, install its own dependencies, and produce the first filled row of the
pre-registered performance envelope before any rare-disease work begins.

## What this notebook executes

| Phase | What runs here | SMOKE profile (default) | FULL profile (Lambda, `SMOKE=False`) |
|---|---|---|---|
| Setup | pinned dependencies installed from PyPI | ✅ | ✅ |
| E1 | frozen snapshot: LINCS L1000 metadata, STRING v12 PPI, 3 GEO bulk series, platform annotations, ChEMBL indications → `data/SNAPSHOT.md` with SHA-256 | ✅ | ✅ |
| E0 | 3-disease calibration panel (Alzheimer's `GSE5281`, rheumatoid arthritis `GSE55235`, pulmonary fibrosis `GSE24206`), disease signatures, connectivity + linear scoring | ✅ real GEO data + synthetic known-answer unit tests of every scorer | ✅ real LINCS Level-5 scoring over all approved perturbagens |
| E4 core | SRCWM latent-dynamics rescue score (cell-type-mixed action dynamics, `min_k` shared objective, `T`-step rollout, `η` scale-matching) | ✅ synthetic world with planted shared-rescue signal, LOGO + bootstrap pipeline | ✅ real roll-outs over LINCS instances, LOGO folds |
| E3 | lifting operator Λ: bulk → per-cell-type via CZ CELLxGENE Census reference | ✅ real brain/lung/blood slices | ✅ larger slices |
| E9 | negative controls: label scrambling, sign-flip, cell-type-affinity scramble | ✅ real data | ✅ |
| E10 | reproducibility artifacts: config hash, manifest, results JSON, `SNAPSHOT.md` | ✅ | ✅ |

## Honest scope statement

* This notebook executes **E0 toolchain validation plus smoke-scale E3/E4 only**.
  The pre-registered full protocol (E2: 5–8 diseases; E4: LOGO × 5–8 folds incl.
  M-CPA / M-GEARS / M-SCGPT / M-STATE arms; E5: 10,000-resample paired bootstrap;
  E6–E8) needs the compute budget of E11 and is a follow-up run, not something
  this file pretends to deliver.
* Deep perturbation models are **not installed or executed** here: honest zero-shot
  adapters require model-specific checkpoint artefacts and considerable install
  steps; faking them in a notebook would violate the design's rules. They are
  described (not executed) in the final section.
* Every number printed comes from data actually downloaded at runtime. Nothing is
  hard-coded. Where a pre-registered gate fails, the notebook prints that it failed.
* No therapeutic claim. Hypothesis-generating scope only (design principle 5).

## Runtime / resources

| Profile | Disk | RAM | Wall time |
|---|---|---|---|
| SMOKE (default) | < 1 GB | < 8 GB | ~15–40 min, CPU only |
| FULL | ~60 GB (LINCS L5 GCTX alone is 21 GB) | 30+ GB | hours |

On the Lambda server set `SMOKE = False` in the configuration cell below.
