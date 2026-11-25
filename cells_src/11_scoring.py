# Connectivity scoring — rank-based KS enrichment (original CMap formulation),
# plus the M-LIN and M-RIDGE linear baselines from EXPERIMENTAL_DESIGN.md E4.
#
# Convention (used consistently everywhere below): larger = better rescue.

def gene_sets(delta, top_k=TOP_K):
    d = delta.sort_values()
    return list(d.index[-top_k:]), list(d.index[:top_k])

def ks_es(ranked_genes, target_set):
    n = len(ranked_genes)
    hit_mask = np.isin(ranked_genes, list(target_set))
    nh = int(hit_mask.sum())
    if nh == 0 or nh == n:
        return 0.0
    step = np.where(hit_mask, 1.0 / nh, -1.0 / (n - nh))
    run = np.cumsum(step)
    k = int(np.argmax(np.abs(run)))
    return float(run[k])

def connectivity_rescue(drug_z, delta, top_k=TOP_K):
    up, dn = gene_sets(delta, top_k)
    order = np.argsort(-drug_z.values, kind="stable")
    ranked = list(drug_z.index[order])
    es_up = ks_es(ranked, up)
    es_dn = ks_es(ranked, dn)
    return float(0.5 * (es_dn - es_up))

def mlin_rescue(drug_z, delta, top_k=TOP_K):
    up, dn = gene_sets(delta, top_k)
    return float(drug_z[dn].mean() - drug_z[up].mean())

def program_cosine(delta, drug_z, background, top_k_pca=PCA_PROGRAMS, seed=SEED):
    V = np.linalg.svd(background - background.mean(axis=1, keepdims=True),
                      full_matrices=False)[2][:top_k_pca]
    dv = (delta - delta.mean()).values @ V.T
    zv = (drug_z - drug_z.mean()).reindex(delta.index).values @ V.T
    denom = (np.linalg.norm(dv) * np.linalg.norm(zv) + 1e-12)
    return float(-(dv @ zv) / denom)

def _test_scores():
    genes = [f"G{i}" for i in range(4000)]
    delta = pd.Series(0.0, index=genes)
    delta.iloc[:1000] = rng.normal(1.0, 0.01, 1000)
    delta.iloc[1000:2000] = rng.normal(-1.0, 0.01, 1000)
    reversing = delta.copy() * -1 + rng.normal(0, 0.01, 4000)
    aggravating = delta.copy() * 1 + rng.normal(0, 0.01, 4000)
    noise = pd.Series(rng.normal(0, 1, 4000), index=genes)
    up, dn = gene_sets(delta)
    for label, z in [("reversing", reversing), ("aggravating", aggravating), ("null", noise)]:
        cs = connectivity_rescue(z, delta)
        ml = mlin_rescue(z, delta)
        print(f"  {label:11s} CS={cs:+.3f}  M-LIN={ml:+.3f}")
    assert connectivity_rescue(reversing, delta) > 0.5
    assert connectivity_rescue(aggravating, delta) < -0.5
    assert abs(connectivity_rescue(noise, delta)) < 0.15
    assert mlin_rescue(reversing, delta) > mlin_rescue(aggravating, delta)
    print("scorer unit tests passed (reversal/aggravation/null separation)")

_test_scores()
