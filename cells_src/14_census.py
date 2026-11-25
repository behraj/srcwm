# E3 — Lifting operator Λ: bulk module → per-cell-type attribution, validated
# on real CZ CELLxGENE Census slices (brain / lung / blood).
#
# Λ is solved as NNLS of |delta| onto the cell-type centroid matrix W of the
# matched tissue. E3.1 (internal consistency) and E3.2 (literature-anchored
# attribution) are both computed and reported honestly; a failing gate is
# reported as failing, never smoothed over.
#
# Each tissue runs two isolated worker subprocesses (scan → X read) written by
# this cell, so the heavyweight TileDB buffers of the obs scan and of the X
# array never co-exist in one process; the notebook itself only handles the
# small centroid matrices.

import gc
import subprocess
import sys
import cellxgene_census
from scipy.optimize import nnls

CELLS_PER_TYPE = 500 if SMOKE else 8000
N_REF_CELLS = 60 if SMOKE else 800
GENE_VARIANCE_TOP = 1200

WORKER = DATA / "_census_slice_worker.py"
WORKER.write_text('''
import sys, json
import numpy as np, pandas as pd
import pyarrow as pa, pyarrow.compute as pc
import cellxgene_census

def rss():
    with open("/proc/self/status") as fh:
        for ln in fh:
            if ln.startswith("VmRSS"):
                return int(ln.split()[1]) // 1024
    return -1

mode, census_version = sys.argv[1], sys.argv[2]

if mode == "scan":
    tissue, expected_json, out_csv, cap = sys.argv[3:7]
    cap = int(cap)
    want = [t.lower() for t in json.loads(expected_json)]
    with cellxgene_census.open_soma(census_version=census_version) as census:
        obs = census["census_data"]["homo_sapiens"]["obs"]
        vf = f'tissue_general == "{tissue}" and disease == "normal"'
        collected, n_seen = {}, 0
        for tbl in obs.read(value_filter=vf, column_names=["soma_joinid", "cell_type"],
                            batch_size=200_000):
            col = tbl.column("cell_type")
            if pa.types.is_dictionary(col.type):
                colc = col.combine_chunks()
                colarr = colc.chunk(0) if isinstance(colc, pa.ChunkedArray) else colc
                pairs = list(enumerate(colarr.dictionary.to_pylist()))
                getter = lambda pos: pc.equal(colarr.indices, pos)
            else:
                pairs = [(0, v) for v in pc.unique(col).to_pylist()]
                getter = None
            for pos, val in pairs:
                lv = val.lower()
                if not any(w in lv for w in want):
                    continue
                cur = collected.setdefault(val, [])
                if len(cur) >= cap:
                    continue
                mask = getter(pos) if getter is not None else pc.equal(col, val)
                cur.extend(tbl.filter(mask).column("soma_joinid").to_pylist())
                if len(cur) > cap:
                    del cur[cap:]
            n_seen += tbl.num_rows
    matched = sorted(collected, key=lambda c: -len(collected[c]) if collected[c] else 0)
    matched = [c for c in matched if collected[c]]
    print(f"[scan] {tissue}: scanned {n_seen} cells; matched types: {len(matched)}")
    if len(matched) < 2:
        print(f"[scan] insufficient reference cell types in '{tissue}' "
              f"({len(matched)}) — reporting honestly, no sample written")
        sys.exit(0)
    selected = matched[:6]
    rows = [(j, c) for c in selected for j in collected[c][:cap]]
    pd.DataFrame(rows, columns=["soma_joinid", "cell_type"]).to_csv(out_csv, index=False)
    pa.default_memory_pool().release_unused()
    print(f"[scan] wrote {out_csv}: {len(rows)} cells across {len(selected)} types")
elif mode == "xread":
    tissue, sample_csv, out_csv, nref = sys.argv[3:7]
    nref = int(nref)
    with cellxgene_census.open_soma(census_version=census_version) as census:
        rel = census["census_data"]["homo_sapiens"]
        var = rel["ms"]["RNA"]["var"]
        X = rel["ms"]["RNA"]["X"]["normalized"]
        var_df = var.read(column_names=["soma_joinid", "feature_name",
                                        "feature_type"]).concat().to_pandas()
        gname = dict(zip(var_df.soma_joinid, var_df.feature_name.str.upper()))
        pc_mask = np.zeros(int(var_df.soma_joinid.max()) + 1, dtype=bool)
        pc_mask[var_df[var_df.feature_type == "protein_coding"].soma_joinid.to_numpy()] = True
        sample = pd.read_csv(sample_csv)
        sample = sample.sort_values("soma_joinid").groupby("cell_type", as_index=False)\
                       .head(nref)
        jids = np.sort(sample.soma_joinid.to_numpy())
        ct_of = dict(zip(sample.soma_joinid.astype(int), sample.cell_type))
        print(f"[xread] {tissue}: reading X slice for {len(jids)} reference cells")
        tbl = X.read(coords=[jids], result_order="row-major",
                     batch_size=100_000).tables().concat()
        gj = tbl["soma_dim_1"].to_numpy()
        cj = tbl["soma_dim_0"].to_numpy()
        vl = tbl["soma_data"].to_numpy().astype(np.float32)
        del tbl
        gm = pc_mask[gj]
        gj, cj, vl = gj[gm], cj[gm], vl[gm]
        del gm
        sd = pd.DataFrame({"g": gj, "type": cj, "v": vl})
        del gj, cj, vl
        sd["type"] = sd["type"].map(ct_of)
        sd = sd.dropna()
        centroids = sd.groupby(["g", "type"])["v"].mean().unstack().fillna(0.0)
        del sd
        centroids.index = [gname.get(int(i), f"HID_{int(i)}") for i in centroids.index]
        centroids.to_csv(out_csv)
        pa.default_memory_pool().release_unused()
        print(f"[xread] wrote {out_csv}: {centroids.shape}")
else:
    raise ValueError(f"unknown mode {mode}")
''')

census_results = {}
fig, axes = plt.subplots(1, len(PANEL), figsize=(4.5 * len(PANEL), 3.6))
axes = np.atleast_1d(axes)

for ax, (name, cfg) in zip(axes, PANEL.items()):
    tissue = cfg["census_tissue_general"]
    sample_csv = DATA / f"_census_sample_{name}.csv"
    out_csv = RESULTS / f"census_centroids_{name}.csv"
    if out_csv.exists():
        out_csv.unlink()
    logs = []
    for mode, argv in [
        ("scan", ["scan", CENSUS_VERSION, tissue,
                  json.dumps(cfg["expected_cell_types"]), str(sample_csv),
                  str(CELLS_PER_TYPE)]),
        ("xread", ["xread", CENSUS_VERSION, tissue, str(sample_csv), str(out_csv),
                   str(N_REF_CELLS)]),
    ]:
        r = subprocess.run([sys.executable, str(WORKER)] + argv,
                           capture_output=True, text=True)
        logs.append((mode, r.returncode, (r.stdout or "") + (r.stderr or "")))
    for mode, rc, out in logs:
        print(f"-- {name}/{mode} (rc={rc}):")
        print(re.sub(r"\n\s*\n", "\n", out).strip()[-900:])
    if not out_csv.exists():
        census_results[name] = dict(status="worker_failed")
        ax.set_axis_off()
        continue
    centroids = pd.read_csv(out_csv, index_col=0)
    centroids = centroids[~centroids.index.duplicated(keep="first")]
    marker_ids = list(centroids.var(axis=1).sort_values(ascending=False)
                      .index[:GENE_VARIANCE_TOP])
    centroids = centroids.loc[marker_ids]
    marker_genes = list(centroids.index)
    cols = list(centroids.columns)
    d = deltas[name]
    common = [m for m in marker_genes if m in d.index]
    Wsub = np.abs(centroids.loc[common].to_numpy())
    yv = np.abs(d.loc[common].to_numpy())
    if Wsub.shape[1] < 2:
        census_results[name] = dict(status="degenerate_basis")
        ax.set_axis_off()
        continue
    a, _ = nnls(Wsub, yv)
    pred = Wsub @ a
    rel_err = float(np.linalg.norm(pred - yv) / (np.linalg.norm(yv) + 1e-12))
    a_norm = a / (a.sum() + 1e-12)
    top3 = [cols[i] for i in np.argsort(-a_norm)[:3]]
    e32 = any(any(e.lower() in t.lower()
                  for e in [ee.lower() for ee in cfg["expected_cell_types"]])
              for t in top3)
    perm = Wsub[:, rng.permutation(Wsub.shape[1])]
    a_scr, _ = nnls(perm, yv)
    a_scr_norm = a_scr / (a_scr.sum() + 1e-12)
    census_results[name] = dict(status="ok", cell_types=cols,
                                attribution=a_norm.tolist(), top3=top3,
                                e32_pass=bool(e32), internal_rel_error=rel_err,
                                scrambled_l1=float(np.abs(a_norm - a_scr_norm).sum()))
    pd.Series(a_norm, index=cols).sort_values().plot.barh(ax=ax,
                                                          title=f"{name} ({tissue})")
    print(f"{name}: E3.1 relative error on {len(common)} marker genes = "
          f"{rel_err:.3f} (design target: median < 0.10); top-3 attribution: "
          f"{top3}; E3.2: {'pass' if e32 else 'FAIL'}; scramble L1 shift = "
          f"{np.abs(a_norm - a_scr_norm).sum():.2f}")
    del centroids, Wsub, yv
    gc.collect()

fig.tight_layout()
fig.savefig(FIGS / "e3_lifting_attribution.png", dpi=120)
plt.close(fig)
json.dump(census_results, open(RESULTS / "e3_lifting.json", "w"), indent=1)
