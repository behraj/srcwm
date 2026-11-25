# SRCWM: Experimental Design

**Project.** Shared-Rescue Cell World Model (SRCWM): cell-type-resolved, multi-disease drug repurposing via action-conditioned latent dynamics with honest zero-shot evaluation.

**Purpose of this document.** Specify every experiment: data, splits, models, baselines, metrics, ablations, compute, and decision gates, so that results are pre-registered before any run and the "when-it-beats-linear" claim is decided by protocol rather than post-hoc selection.

---

## 0. Design Principles

1. **Pre-registration.** All splits, hyperparameter grids, and primary metrics are fixed in configuration files before the first run. Post-hoc metric additions are labelled as exploratory in every table.
2. **Linear baselines are primary.** Every headline question is answered as *SRCWM vs. linear mean-shift/ridge*, not SRCWM vs. other deep models. Beating weak baselines only is a negative result and is reported as such.
3. **Zero-shot by construction.** Any statistic involving disease \(d\) or cell type \(c\) must not use $d$s labels or $c$s calibration data anywhere in training, lifting calibration, or threshold fitting.
4. **Fail informatively.** Negative and mixed outcomes feed the regime-analysis contribution (Experiment E6); nothing is discarded.
5. **Hypothesis-generating scope.** No claim of therapeutic efficacy; every candidate is a testable hypothesis with a mechanistic narrative.

---

## E0. Calibration Study (Milestone 0)

**Goal.** Validate the toolchain on a *known-answer* setting and estimate the linear/deep performance envelope before touching rare diseases.

- **Setting.** 3 common-model diseases with well-curated LINCS reversal results (e.g., inflammatory (rheumatoid arthritis), neurodegenerative (Alzheimer's), oncology-adjacent (fibrosis/uveitis)) with known repurposing outcomes in the literature.
- **Procedure.** Reproduce: (a) standard CMap connectivity scoring; (b) linear mean-shift ranking; (c) one deep perturbation model's zero-shot score. Compare against known indications.
- **Deliverables.**
  - Toolchain sanity: identical gene-space mapping, identical signature processing across sources.
  - A first filled "performance envelope" row: where does each method rank known indications?
  - Scale-matching constant \(\eta\) (latent step size) and calibration thresholds \((r_0, w_c)\).
- **Gate.** Proceed only if linear reproduction reaches published-level recall of known indications. If not, fix the data pipeline first; do not tune the model.

---

## E1. Data Infrastructure

### Datasets


| Component                 | Source                                        | Onto/ID harmonization                                        | Minimal version                    |
| ------------------------- | --------------------------------------------- | ------------------------------------------------------------ | ---------------------------------- |
| Disease bulk signatures   | GEO                                           | Entrez → Ensembl → gene symbol; reanalyse raw where possible | 5 diseases, ≥ 8 case/control each  |
| Drug signatures           | LINCS L1000 / CMap                            | CMap gene symbols; replicate QC (mod-z)                      | approved drugs, ≥ 3 replicates     |
| Single-cell reference     | public atlases (healthy tissue-matched)       | same gene space; cell-type labels harmonized                 | one reference per disease tissue   |
| Cell-type marker matrices | atlas + literature                            | curated \(W_c\)                                              | per disease tissue                 |
| PPI                       | STRING (high confidence)                      | protein → gene symbols                                       | one release, frozen                |
| Targets / indications     | ChEMBL, DrugBank, Open Targets                | ChEMBL-ID drug keys                                          | frozen snapshot, dated in repo     |
| Perturb-seq (optional)    | public datasets (e.g., Replogle/Norman-style) | gene space mapping                                           | one dataset for world-model sanity |


### Frozen snapshots

Every database is fetched once into an immutable data lake with a `SNAPSHOT.md` (date, version, checksums). All experiments cite one snapshot.

---

## E2. Disease Panel Selection (Milestone 1, Part A)

**Selection criteria (stated before selection):**

1. GEO bulk data with ≥ 8 cases/controls and public raw availability.
2. Tissue covered by a single-cell atlas with ≥ 3 resolvable cell types (enables lifting).
3. Clinically heterogeneous categories: include at least one neurodegenerative, one metabolic/renal, one immunologic, one rare developmental, one hematologic.
4. At least 2 diseases with documented shared pathways for closed-loop validation, and at least 2 with no literature-documented shared therapy (genuine discovery space).
5. No disease appears in any training split of the world model's pretraining corpus (checked against LINCS disease-tagged perturbations).

**Shortlist procedure.** Score candidate diseases against criteria 1–5, take the top 5–8. Final panel recorded with a written justification per disease *before* any model run.

---

## E3. Lifting Operator Validation (Milestone 1, Part B)

**Goal.** Show \(\Lambda(\delta_k)\) (bulk → per-cell-type) is faithful before any roll-outs depend on it.

- **E3.1 Internal consistency.** Re-aggregate lifted cell-type deltas and confirm recovery of \(\delta_k\) (relative error target: median &lt; 10% on disease genes).
- **E3.2 External consistency.** For each disease with disease-relevant cell types known independently from literature, the lifting mass attribution must rank those cell types within the top 3 (per-diseased basis; report successes and failures).
- **E3.3 Ablation.** Compare lifting against naive tissue-level injection (\(\eta\) fixed, no cell-type prior) on downstream rescue scoring for one pair of diseases; the lift must agree or exceed.
- **Gate.** If E3.2 fails on &gt; 2 diseases, revisit basis matrices or cell-type annotation before proceeding.

---

## E4. Shared-Rescue Roll-Outs (Milestone 2)

### Models under test


| ID      | Model                          | Type            | Frozen?                |
| ------- | ------------------------------ | --------------- | ---------------------- |
| M-LIN   | linear mean-shift / CMap score | baseline        | n/a                    |
| M-RIDGE | ridge on gene-program PCs      | baseline        | n/a                    |
| M-CPA   | CPA / chemCPA                  | deep            | yes, + disease adapter |
| M-GEARS | GEARS                          | deep            | yes, + disease adapter |
| M-SCGPT | scGPT + action head            | deep/foundation | yes, + adapter         |
| M-STATE | Arc Institute State            | deep/foundation | yes, + adapter         |


### Protocol

1. **Split by disease:** one disease held fully out per evaluation (leave-one-disease-out, LOGO). Panel of 5 → 5 folds per split; panel of 8 → 8 folds.
2. **Split by cell type:** within each disease, hold out the 2 cell types with the strongest residual pathology from any per-cell-type calibration.
3. **Drug ranking.** For each fold, score all approved drugs with each model's rescue score. No threshold tuning on held-out data.
4. **Multi-disease composition.** SRCWM's panel-level \(\min_k\) score is computed per LOGO fold (i.e., excluding the held-out disease from the objective), so generalization to the held-out disease tests *transfer of the shared-rescue hypothesis*, not memorization.

### Primary metrics (fixed pre-run)


| Metric                  | Definition                                                                                       | Justification                                                     |
| ----------------------- | ------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------- |
| P@25                    | fraction of top-25 ranked drugs with known indication or literature support for held-out disease | early-precision (a scarce-resource ranking criterion)             |
| NDB                     | normalized discount-benefit of known indications vs. geometric controls                          | power tradeoff                                                    |
| Dose-consistency        | Spearman \(\rho\)(score, dose) over L1000 dose curves of top-50                                  | mechanistic sanity (signature reversal should be dose-consistent) |
| Program-rescue coverage | fraction of disease-active programs positively rescued by top-25                                 | explanation quality                                               |
| Safety-margin diversity | number of distinct target narratives among top-25                                                | signals multiple mechanisms rather than one metric artifact       |


### Secondary metrics (reported but labelled secondary)

AUROC over all ranked drugs, mean reciprocal rank, computational cost per scoring run.

---

## E5. Honest Baseline Comparison (Milestone 2 cont.)

**Headline question.** Does SRCWM beat M-LIN/M-RIDGE on P@25 and NDB under LOGO splits, and by how much?

- Paired bootstrap over drugs (10,000 resamples), Holm-corrected across the 5 LOGO folds.
- **Pre-registered decision rule:** if SRCWM does not beat M-LIN on ≥ 3 of 5 folds at \(p < 0.05\) (Holm), the primary conclusion is negative and the paper's emphasis shifts to the regime analysis (E6). This rule is binding; no post-hoc switch to secondary metrics as headline.
- Additionally report a **chance-calibrated control**: reshuffle drug identity within dose buckets (placebo ranking) to confirm positive signal's magnitude.

---

## E6. When-It-Beats-Linear Regime Analysis (the methodological contribution)

Stratify every LOGO result along three axes and analyze where SRCWM's advantage concentrates:

1. **Disease novelty:** number of diseases in the panel most similar (pathway overlap with STRING/Jaccard) to the held-out disease.
2. **Cell-type novelty:** fraction of held-out disease's affected cell types absent from calibration.
3. **Dose density:** number of available LINCS replicate doses for each candidate drug.

Deliverable: a 3-panel "advantage map" with confidence bands, plus a written statement of regimes where the world model adds value vs. regimes where it does not. This section is written with the same care as a positive result and can stand alone as the paper's contribution if effects are modest.

---

## E7. Counterfactual Attribution and Mechanistic Validation (Milestone 3)

**Goal.** Show explanations are anchored, not decorative.

### E7.1 Attribution faithfulness

- Perturb latent roll-out inputs (drop top-attributed programs, re-score). Predicted rescue should degrade proportionally to the removed attribution mass (deletion test).
- Report Spearman correlation between attribution mass and deletion-induced score change across top-50 candidates (target \(\rho > 0.4\)).

### E7.2 Target anchoring

- For each top-25 drug per held-out disease: fraction of positively-rescued programs reachable from the drug's primary targets within 2 STRING hops with FDR-corrected enrichment \(< 0.05\).
- Compare against a target-shuffled null (same drug set, permuted targets, 1000 draws).

### E7.3 Retrospective literature validation

- Cross-reference each top-25 mechanism narrative against ChEMBL, Open Targets, and PubMed (structured search template released with the paper).
- Report three counts: confirmed mechanisms, contradicted mechanisms, and uncharacterized; the last group defines the experimental-validation queue.

### E7.4 Docking (optional, flagged as mechanistic corroboration only)

- For the top 5 shared-rescue candidates with well-resolved protein structures, perform molecular docking of the drug against the primary target; report affinity ranks relative to decoy ligands. Included only as corroboration, never as evidence of efficacy.

---

## E8. Ablation Grid


| Ablation                              | Description                            | Expected signal                                                              |
| ------------------------------------- | -------------------------------------- | ---------------------------------------------------------------------------- |
| A1 lifting off                        | tissue-level injection, \(\eta\) fixed | degradation concentrated on multi-cell-type diseases                         |
| A2 \(\min_k\) → mean                  | shared objective removed               | top-ranked drugs become single-disease specialists; multi-disease P@25 falls |
| A3 frozen vs. fine-tuned              | unfreeze full model                    | risk of overfitting to calibration diseases; report both                     |
| A4 roll-out depth \(T \in \{1,2,3\}\) | depth sweep                            | T=1 primary; &gt;1 exposes accumulation error                                |
| A5 \(\lambda_2\) sweep                | cell-type prior strength \([0,10]\)    | U-shaped or monotonic sensitivity; report                                    |
| A6 model swap                         | scGPT vs. State vs. CPA                | robustness of the shared-rescue conclusion to world-model choice             |


Run the full grid on the 2 most-documented diseases and the widest grid on the remainder, to keep compute tractable.

---

## E9. Negative Controls


| Control                                         | Purpose                                                                                                             |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| Scrambled tissue labels (case/control shuffled) | should suppress disease signal (all metrics → chance)                                                               |
| Scrambled cell-type affinities                  | isolated the cell-type prior's contribution                                                                         |
| Signature sign-flip                             | should *lower* the score if score truly detects reversal, not signal magnitude alone (magnitude-collapse detection) |
| Plate / batch effects                           | within-LINCS plate shuffling; metrics must not recover under plate-only signal                                      |
| Dose-directionality                             | for drugs with monotone dose responses, score must be monotone in dose                                              |


---

## E10. Reproducibility Contract

- One repository: `srcwm/` with `configs/`, `data/SNAPSHOT.md`, `experiments/`, `results/`.
- Every experiment identified by a hash of (config, snapshot, split seed). Results keyed by that hash; no duplicated manual runs.
- All splits, hyperparameter grids, and primary metrics in version-controlled YAML, committed *before* the corresponding result directory exists (enforce via CI timestamp check).
- Renders of every figure generated by a single `make figures` target in a clean container.

---

## E11. Compute and Timeline


| Item                                          | Estimate                                                                                  |
| --------------------------------------------- | ----------------------------------------------------------------------------------------- |
| Frozen model inference (all folds, all drugs) | CPU-bound: \~1–2 GPU-days total on a single 24-GB GPU; embarrassingly parallel over drugs |
| Adapter training                              | &lt; 4 GPU-hours per LOGO fold                                                            |
| Lifting operator                              | CPU, minutes per disease                                                                  |
| Docking (optional E7.4)                       | CPU cluster, \~1–2 days for top-5                                                         |



| Phase                                    | Duration    | Gate                                    |
| ---------------------------------------- | ----------- | --------------------------------------- |
| P0: toolchain + calibration (E0, E1)     | weeks 1–4   | E0 passes                               |
| P1: panel + lifting (E2, E3)             | weeks 5–10  | E3.2 ≥ 1-holdout success rate           |
| P2: roll-outs + baselines (E4, E5)       | weeks 11–18 | E5 pre-registered decision rule applied |
| P3: attribution + ablations (E7, E8, E9) | weeks 19–26 | attribution correlation target met      |
| P4: regime analysis + writing (E6, E10)  | weeks 27–32 | full reproducibility contract met       |


*Risk buffer of 4–6 weeks assumed for data-quality surprises (inconsistent GEO processing, atlas co-coverage gaps).*

---

## 12. Risk Register


| Risk                                                                              | Likelihood  | Mitigation                                                                                                  |
| --------------------------------------------------------------------------------- | ----------- | ----------------------------------------------------------------------------------------------------------- |
| Deep model fails to beat linear on any fold                                       | medium      | pre-registered pivot to E6 regime analysis as headline                                                      |
| Lifting operator attributes mass to wrong cell types                              | medium      | E3.2 gate; fall back to marker-gene-constrained basis                                                       |
| Zero-shot contamination via foundation-model pretraining                          | low/unknown | inspect per-disease and per-drug overlap in the world model's training card; document any suspected leakage |
| Too few approved drugs with LINCS coverage to make multi-disease stats meaningful | low         | widen to investigational compounds; report selection sensitivity                                            |
| Reproducibility of 2026-dated preprints (API/data changes)                        | medium      | frozen snapshots; walled garden of data versions; document retrieval dates                                  |


---

## 13. Success Criteria (pre-registered)

The project succeeds — regardless of whether SRCWM does or does not outperform linear baselines — if it delivers:

1. A repeatable, honestly split, multi-disease rescue benchmark with released code and data snapshots.
2. A quantified "when-it-beats-linear" regime map for cell world models in repurposing.
3. At least one concrete, mechanistically explained, shared-rescue candidate prioritized for experimental validation, with a residual evidence gap stated plainly.
4. Reproducibility artifacts sufficient for a third party to regenerate every number.
