# E10 — reproducibility artifacts.

import platform, sys, scipy

RUN_INFO = dict(config_hash=CONFIG_HASH, created_utc=CONFIG["created_utc"],
                python=sys.version.split()[0], platform=platform.platform(),
                numpy=np.__version__, pandas=pd.__version__,
                scipy=scipy.__version__, census_version=CENSUS_VERSION,
                seed=SEED, profile="SMOKE" if SMOKE else "FULL")
catalog = {
    "results/e0_disease_deltas.json": "disease signatures (real GEO ingestion)",
    "results/e4_synthetic_known_answer.json": "synthetic world LOGO P@25 / NDB",
    "results/e4_synthetic_bootstrap.json": "paired bootstrap + Holm (synthetic)",
    "results/e3_lifting.json": "Census lifting E3.1 / E3.2 + scramble control",
    "results/e9_label_scramble.json": "E9 label-scramble on real series",
    "results/e0_run_summary.json": "profile summary / gate status",
    "figures/e3_lifting_attribution.png": "per-disease cell-type attribution",
}
present = {}
for rel, desc in catalog.items():
    ok = Path(rel).exists()
    present[rel] = dict(written=ok, description=desc,
                        sha256=sha256_file(rel) if ok else None)
census_slice_catalog = {}
for rel, desc in {
    f"data/_census_sample_{name}.csv":
        f"Census slice selection, {name} ({PANEL[name]['census_tissue_general']})",
    f"results/census_centroids_{name}.csv":
        f"Census cell-type centroids, {name} ({PANEL[name]['census_tissue_general']})",
}.items():
    ok = Path(rel).exists()
    census_slice_catalog[rel] = dict(
        written=ok, description=desc, census_version=CENSUS_VERSION,
        bytes=Path(rel).stat().st_size if ok else None,
        sha256=sha256_file(rel) if ok else None)
json.dump(dict(run=RUN_INFO, artifacts=present,
               census_slices=census_slice_catalog, manifest=MANIFEST),
          open(RESULTS / "e0_reproducibility_bundle.json", "w"), indent=1)
print("results catalog:")
for rel, info in present.items():
    print(f"  {'OK ' if info['written'] else 'MISS'} {rel:45s} {info['description']}")
print("census slice hashes:")
for rel, info in census_slice_catalog.items():
    print(f"  {'OK ' if info['written'] else 'MISS'} {rel:45s} "
          f"{info['sha256'][:16] if info['sha256'] else '-'}")
