from scipy.optimize import nnls

G, C = 978, 4
K_DISEASES = ["D1", "D2", "D3"]
genes = [f"G{i}" for i in range(G)]

W = np.zeros((G, C))
for c in range(C):
    support = rng.choice(G, size=120, replace=False)
    W[support, c] = rng.gamma(2.0, 1.0, 120)
W /= np.linalg.norm(W, axis=0, keepdims=True) + 1e-12

adj = rng.uniform(0, 1, (C, C)) + 0.3 * np.eye(C)
M = adj / adj.sum(axis=1, keepdims=True)

lat_deltas = {}
for k, idx in zip(K_DISEASES, [[0, 1], [1, 2], [2, 3]]):
    a = np.zeros(C)
    for j, wgt in zip(idx, (0.8, 0.5)):
        a[j] = wgt
    lat_deltas[k] = a

delta_genes = {k: pd.Series(W @ a, index=genes) for k, a in lat_deltas.items()}

n_drugs = 120
V = np.zeros((C, n_drugs))
drug_names = [f"DRUG_{i}" for i in range(n_drugs)]
shared = [f"DRUG_{i}" for i in range(3)]
specialists = [f"DRUG_{i}" for i in range(3, 33)]

for j in range(3):
    V[:, j] = -np.mean(np.stack([lat_deltas[k] for k in K_DISEASES]), axis=0) * 1.5 \
              + rng.normal(0, 0.02, C)
for j in range(3, 33):
    c = rng.integers(0, C)
    dists = [abs(lat_deltas[kk][c]) for kk in K_DISEASES]
    k = K_DISEASES[int(np.argmax(dists))]
    V[:, j] = -lat_deltas[k] * 1.2 + rng.normal(0, 0.05, C)
for j in range(33, n_drugs):
    V[:, j] = rng.normal(0, 1.0, C) / np.sqrt(C)

lifted = {k: a / (np.linalg.norm(a) + 1e-12) for k, a in lat_deltas.items()}
ETA = 0.5

def rollout(v, eta=ETA, T=1):
    x = v.copy()
    for _ in range(T):
        x = (1 - eta) * x + eta * (M @ x)
    return x

def srcwm_rescue(drug_v, lifted_deltas, T=1):
    v = rollout(drug_v, ETA, T)
    out = {}
    for k, ld in lifted_deltas.items():
        out[k] = float(-(v @ ld) / (np.linalg.norm(v) * np.linalg.norm(ld) + 1e-12))
    return out

def drug_signature(name_i, noise=0.01):
    z = W @ V[:, name_i] + rng.normal(0, noise, G)
    return pd.Series((z - z.mean()) / (z.std() + 1e-12), index=genes)

src_scores, lin_scores = {}, {}
for i, name in enumerate(drug_names):
    z = drug_signature(i)
    v = W.T @ (z.values - z.values.mean())
    src_scores[name] = srcwm_rescue(v, lifted)
    lin_scores[name] = {k: mlin_rescue(z, delta_genes[k]) for k in K_DISEASES}

def rank_by(scores_per_drug, train, agg="min"):
    key = (lambda d: -min(scores_per_drug[d][k] for k in train)) if agg == "min" \
          else (lambda d: -np.mean([scores_per_drug[d][k] for k in train]))
    return sorted(drug_names, key=key)

def dcg(ranking, truth, k=25):
    return sum(1.0 / np.log2(i + 2) for i, d in enumerate(ranking[:k]) if d in truth)

def p_at_k(ranking, truth, k=25):
    return len(set(ranking[:k]) & set(truth)) / k

def ndb(ranking, truth, k=25):
    hd = sum(1.0 / np.log2(i + 2) for i in range(k))
    rand = hd * len(truth) / n_drugs
    ideal = sum(1.0 / np.log2(i + 2) for i in range(min(k, len(truth))))
    denom = ideal - rand
    return 0.0 if denom <= 0 else (dcg(ranking, truth, k) - rand) / denom

folds = {}
for held in K_DISEASES:
    train = [k for k in K_DISEASES if k != held]
    r_src = rank_by(src_scores, train)
    r_lin = rank_by(lin_scores, train)
    truth = set(shared)
    r_mean = rank_by(src_scores, train, agg="mean")
    folds[held] = dict(P25_srcwm=p_at_k(r_src, truth), P25_mlin=p_at_k(r_lin, truth),
                       NDB_srcwm=ndb(r_src, truth), NDB_mlin=ndb(r_lin, truth),
                       shared_in_top5=sum(s in r_src[:5] for s in shared),
                       specialist_top1_mean_aggregate=r_mean[0])
folds_report = pd.DataFrame(folds).T
print(folds_report.drop(columns=["specialist_top1_mean_aggregate"]))
print("shared-rescue drugs in top-5 of min_k ranking per fold:",
      {k: f["shared_in_top5"] for k, f in folds.items()})

for held, f in folds.items():
    assert f["P25_srcwm"] >= 2 / 25, f"known-answer gate failed on {held}"
print("E0 known-answer gate: planted shared-rescue drugs recovered by SRCWM min_k on all folds")

def holm(pvals):
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order[::-1]):
        running = max(running, (m - rank) * p[idx])
        adj[idx] = min(running, 1.0)
    return adj

boot_rows = []
rho = {}
for held in K_DISEASES:
    train = [k for k in K_DISEASES if k != held]
    r_src = rank_by(src_scores, train)
    r_lin = rank_by(lin_scores, train)
    truth = set(shared)
    src_hits = np.array([1 if d in truth else 0 for d in r_src])
    lin_hits = np.array([1 if d in truth else 0 for d in r_lin])
    top = np.concatenate([np.ones(25), np.zeros(n_drugs - 25)])
    diffs = []
    for b in range(BOOTSTRAP_B):
        idx = rng.integers(0, n_drugs, n_drugs)
        s_hits = src_hits[idx] * top
        l_hits = lin_hits[idx] * top
        diffs.append((s_hits.sum() / 25) - (l_hits.sum() / 25))
    diffs = np.array(diffs)
    se = diffs.std(ddof=1) / np.sqrt(len(diffs))
    from scipy.stats import norm
    p = 2 * (1 - norm.cdf(abs(diffs.mean()) / (se + 1e-12)))
    boot_rows.append(dict(fold=held, delta_P25=diffs.mean(),
                          ci_lo=np.percentile(diffs, 2.5),
                          ci_hi=np.percentile(diffs, 97.5), p=p))
    rho[held] = diffs

pv = holm(np.array([r["p"] for r in boot_rows]))
for r, p_h in zip(boot_rows, pv):
    r["p_holm"] = p_h
boot_df = pd.DataFrame(boot_rows)
print(boot_df)

W_scr = W[rng.permutation(G)]
a_scr, _ = nnls(W_scr, np.abs(delta_genes["D1"].values))
ld_scr = a_scr / (np.linalg.norm(a_scr) + 1e-12)
scrambled_rescue = {name: srcwm_rescue(W_scr.T @ (drug_signature(i).values
                                                  - drug_signature(i).values.mean()),
                                       {"D1": ld_scr})["D1"]
                    for i, name in enumerate(drug_names)}
r_scr = sorted(drug_names, key=lambda d: -scrambled_rescue[d])
P25_scramble = p_at_k(r_scr, set(shared))
print(f"E9 scrambled-cell-type-affinity control on synthetic world: "
      f"P@25 = {P25_scramble:.3f} (unscrambled folds: "
      f"{[folds[k]['P25_srcwm'] for k in folds]})")
assert P25_scramble < 0.08, "scramble control did not collapse — scoring is label-leaky"
print("E9 scrambled-delta control: collapsed as required (score tracks real biology)")

json.dump({k: {kk: float(vv) for kk, vv in f.items() if isinstance(vv, (int, float))}
           for k, f in folds.items()},
          open(RESULTS / "e4_synthetic_known_answer.json", "w"))
json.dump({r["fold"]: {kk: float(vv) for kk, vv in r.items() if kk != "fold"}
           for r in boot_rows},
          open(RESULTS / "e4_synthetic_bootstrap.json", "w"))
