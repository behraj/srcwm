# E9 — negative controls summary (label-scramble on real GEO signatures;
# sign-flip and scramble controls already executed on the scoring and lifting
# machinery above).

rows = []
for name, st in ingested.items():
    g, case, control = st["expr"], st["case"], st["control"]
    n_case, n_ctl = int(case.sum()), int(control.sum())
    labels_all = case | control
    keep_idx = np.flatnonzero(labels_all)
    pm = rng.permutation(int(labels_all.sum())) < n_case
    shuffled_case = np.zeros_like(case)
    shuffled_case[keep_idx] = pm
    h = g.iloc[:, keep_idx].apply(zscore_per_sample)
    case_within = case[keep_idx]
    d_real = h.loc[:, case_within].mean(axis=1) - h.loc[:, ~case_within].mean(axis=1)
    d_scr = h.loc[:, pm].mean(axis=1) - h.loc[:, ~pm].mean(axis=1)
    rows.append(dict(disease=name, case_n=n_case, control_n=n_ctl,
                     real_abs_median=float(d_real.abs().median()),
                     scrambled_abs_median=float(d_scr.abs().median())))
e9 = pd.DataFrame(rows).set_index("disease")
e9["collapse_ratio (scrambled/real)"] = e9.scrambled_abs_median / e9.real_abs_median
print(e9.round(4))
assert (e9["collapse_ratio (scrambled/real)"] < 1).all(), \
    "label scramble must shrink the disease signal (E9 gate)"
print("E9 label-scramble gate: passed on all calibration diseases "
      "(scrambled |Δ| < real |Δ|)")
json.dump(e9.reset_index().to_dict(orient="records"),
          open(RESULTS / "e9_label_scramble.json", "w"))
