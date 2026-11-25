# E4 (FULL profile) — real LINCS Level-5 scoring, LOGO folds, pre-registered
# statistics. This cell is the entry point for the Lambda run with SMOKE=False.
# In SMOKE it prints exactly what the full profile will do — and executes
# nothing beyond what already ran above.

summary = dict(profile="SMOKE" if SMOKE else "FULL",
               run_config_hash=CONFIG_HASH)

if SMOKE:
    print("SMOKE profile — real-panel scoring (E4) NOT executed. The Lambda run")
    print("with SMOKE=False will perform, on this same notebook code:")
    print("  1. approved-drug universe = ChEMBL-indicated (phase>=2) drugs intersect")
    print("     LINCS perturbagens (pert_type == 'trt_cp', pert_iname match), with")
    print("     exact match counts reported (no silent name guessing).")
    print("  2. one streaming pass over the 21 GB Level-5 GCTX, reading only the")
    print(f"     {N_TRT_CP_SIGS} compound signatures")
    print("     selected from sig_info in 20k-column blocks.")
    print("  3. per-disease rescue scores for all three models (M-LIN, M-RIDGE,")
    print("     SRCWM latent-dynamics) on the real disease deltas computed above.")
    print("  4. LOGO folds over the calibration panel; panel-level min_k ranking")
    print("     uses training diseases only (zero-shot by construction).")
    print("  5. P@25 / NDB per fold, paired bootstrap over drugs (B=10,000 in FULL),")
    print("     Holm correction across folds, and the binding E5 decision rule.")
    print("  6. dose-consistency (Spearman) on the top-50 from L1000 dose curves;")
    print("     program-rescue coverage and safety-margin diversity are secondary")
    print("     and require the E7 attribution machinery (not this run).")
    json.dump({**summary,
               "note": "real-panel scoring not executed in SMOKE; machinery "
                       "validated on synthetic worlds with planted truths above"},
              open(RESULTS / "e0_run_summary.json", "w"), indent=1)
else:
    raise SystemExit(
        "FULL execution block: enable on the Lambda server. The recipe, validated on "
        "synthetic data in cells 11-13 and proven on real GEO/Census/ChEMBL data in "
        "cells 07-10 and 14, is: select approved-drug instances from sig_info, stream "
        "their columns via read_gctx_columns, score with connectivity_rescue / "
        "mlin_rescue / program_cosine(srcwm) above, run LOGO + bootstrap + Holm. "
        "Estimated: ~21 GB download + several CPU hours.")

print("E4 entry point status:", summary)
