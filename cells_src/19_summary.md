# Run summary and how to go to FULL on the Lambda server

## What actually executed in this run

| Item | Status |
|---|---|
| Dependency install (cell 1) | executed, versions printed |
| E1 frozen snapshot downloads + LINCS SHA512 verification | executed against live FTP |
| ChEMBL indication snapshot (paged REST) | executed against live API |
| E0 panel ingestion: `GSE5281` (AD), `GSE55235` (RA), `GSE24206` (IPF) | executed, real GEO data |
| Scorer unit tests (CMap connectivity, M-LIN) | executed, asserted |
| GCTX streaming reader | executed on synthetic round-trip (SMOKE) |
| E4 known-answer world: SRCWM vs M-LIN, LOGO + bootstrap + Holm | executed, asserted |
| E3 lifting on real Census slices + E9 scramble | executed against live Census |
| E9 label scramble on real GEO | executed |
| E2 panel shortlist (5–8 diseases), E6–E8, full E4 | **NOT executed** (out of E0 scope) |
| Deep arms M-CPA/M-GEARS/M-SCGPT/M-STATE | **NOT executed** (see stub above) |

## Going from SMOKE to FULL on your Lambda server

1. `git clone https://github.com/behraj/srcwm && cd srcwm`
2. Launch the notebook with an editor/terminal and set, in the config cell:
   `SMOKE = False`
3. Run all cells top-to-bottom, CPU is fine for E0/E1/E3; the real E4 loop is
   CPU-bound and embarrassingly parallel over drug signatures (E11 table).
4. Disk: ~25 GB free for metadata + Level-5 GCTX + slices, plus the 21 GB
   compressed GCTX (~60 GB recommended in for safety while streaming).
5. All outputs land in `results/` and `figures/`; E10 bundle
   `results/e0_reproducibility_bundle.json` carries the config hash, runtime
   versions, and artifact checksums for the reproducibility contract.

## Honesty ledger

- Every number above was computed at runtime from primary-source downloads.
- All pre-registered gates in SMOKE profile either passed (printed above) or
  their failure is printed verbatim — the notebook aborts, it does not skip.
- The FULL-profile E4 loop over all LINCS perturbagens is enabled but its body
  only prints the validated recipe; running it is the P2 milestone (weeks 11–18),
  after the E0 gate and the E3 gates have been reviewed.
