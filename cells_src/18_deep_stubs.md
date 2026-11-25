# Deep-model arms (E4) — specification stub, deliberately **not executed**

The pre-registered E4 protocol includes four deep perturbation-model arms as
frozen models + disease adapters:

| ID | Model | Zero-shot adapter plan | Install footprint |
|---|---|---|---|
| M-CPA | CPA / chemCPA | freeze encoder/decoder, fit disease embedding on calibration panel only | PyTorch + `chemCPA` repo pin |
| M-GEARS | GEARS | freeze gene co-expression model, adapter on disease delta | PyTorch + `gears` (PyG) |
| M-SCGPT | scGPT + action head | frozen foundation weights, linear action head on panel | PyTorch + checkpoint (~100 MB–1 GB) |
| M-STATE | Arc Institute State | frozen foundation model, disease adapter | PyTorch + checkpoint |

**Why this notebook does not run them (honesty):**

1. E0's pre-registered scope is the **linear reproduction gate**: "proceed only if
   linear reproduction reaches published-level recall of known indications". The
   deep arms are only *added after* that gate passes.
2. Each arm needs model-specific checkpoint artefacts (megabyte+ downloads,
   outside public PyPI versioning) and install steps that can break a notebook;
   silently degrading to a stub implementation would fake a headline number and
   violate design principle 2 and the E5 binding decision rule.
3. The E5 comparison is **SRCWM vs. M-LIN** — the deep arms matter for E6 regime
   analysis, not for the Milestone-0 gate.

**TODO (post-E0, on Lambda, in dedicated scripts — not this notebook):**

- [ ] Pin exact commits/tags for chemCPA, GEARS, scGPT, State in `configs/models.yaml`.
- [ ] Download checkpoints once into the frozen data lake, checksum into `SNAPSHOT.md`.
- [ ] Implement per-model `rescue(drug_id, panel) -> {disease: float}` adapters.
- [ ] Re-run E4/E5 with all 6 model-IDs; every fold's P@25/NDB lands in `results/e5_{config_hash}.json`.
