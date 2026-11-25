def parse_series_matrix(path):
    txt = gzip.decompress(Path(path).read_bytes()).decode(errors="ignore")
    sample_ids = [s.strip('"') for s in
                  txt.split("!Sample_geo_accession\t")[1].split("\n")[0].split("\t")[1:]]
    label_rows = []
    for line in re.findall(r'!Sample_characteristics_ch1\t(.*)', txt):
        label_rows.append([v.strip('"') for v in line.split("\t")[1:]])
    t0 = txt.index("!series_matrix_table_begin")
    t1 = txt.index("!series_matrix_table_end")
    body = re.sub(r'^!?series_matrix_table_(begin|end).*$', '', txt[t0:t1], flags=re.M)
    expr = pd.read_csv(io.StringIO(body), sep="\t", index_col=0)
    meta = pd.DataFrame(label_rows).T
    meta.index = sample_ids
    return expr, meta

def parse_annot(path):
    txt = gzip.decompress(Path(path).read_bytes()).decode(errors="ignore")
    i = txt.index("!platform_table_begin")
    j = txt.index("!platform_table_end")
    rows = txt[txt.index("\n", i) + 1:j].strip().split("\n")
    header = rows[0].split("\t")
    df = pd.read_csv(io.StringIO("\n".join(rows[1:])), sep="\t", header=None,
                     names=header, usecols=["ID", "Gene symbol"], dtype=str,
                     engine="c").dropna()
    df["Gene symbol"] = df["Gene symbol"].str.split("///").str[0].str.strip().str.upper()
    return dict(zip(df.ID, df["Gene symbol"]))

def zscore_per_sample(s):
    return (s - s.mean()) / (s.std() + 1e-12)

def sample_labels(meta, sid):
    return " | ".join(str(v) for v in meta.loc[sid].values).lower()

def extract_region(meta, sid, contains):
    vals = [str(v).lower() for v in meta.loc[sid].values]
    for v in vals:
        if contains in v:
            return v.split(":", 1)[-1].strip()
    return "unknown"

def ingest(name, cfg):
    acc = cfg["acc"]
    expr, meta = parse_series_matrix(DATA / "geo" / f"{acc}_series_matrix.txt.gz")
    sym = parse_annot(DATA / "geo" / f"{cfg['gpl']}.annot.gz")
    no_meta = [c for c in expr.columns if c not in meta.index]
    if no_meta:
        print(f"  {acc}: dropping {len(no_meta)} sample column(s) lacking series "
              f"metadata: {no_meta[:5]}")
        expr = expr[[c for c in expr.columns if c in meta.index]]
    p2s = pd.Series([sym.get(p) for p in expr.index], index=expr.index).dropna()
    g = expr.loc[p2s.index]
    g.index = p2s.values
    g = g.groupby(level=0).mean()
    labels = [sample_labels(meta, sid) for sid in g.columns]
    case = np.array([any(c in labels[i] for c in [t.lower() for t in cfg["case_terms"]])
                     for i in range(g.shape[1])])
    control = np.array([any(c in labels[i] for c in [t.lower() for t in cfg["control_terms"]])
                        for i in range(g.shape[1])])
    if cfg.get("stratify") == "region":
        zs = g.apply(zscore_per_sample)
        region_of = [extract_region(meta, sid, cfg["region_line_contains"])
                     for sid in g.columns]
        per_region = {}
        for reg in sorted(set(region_of)):
            cols = [c for c, r in zip(g.columns, region_of) if r == reg]
            cm = np.array([any(c in sample_labels(meta, cc).lower()
                               for c in [t.lower() for t in cfg["case_terms"]])
                           for cc in cols])
            km = np.array([any(c in sample_labels(meta, cc).lower()
                               for c in [t.lower() for t in cfg["control_terms"]])
                           for cc in cols])
            if cm.sum() >= 1 and km.sum() >= 1:
                zs_reg = zs[cols]
                per_region[reg] = zs_reg.T[cm].mean(axis=0) - zs_reg.T[km].mean(axis=0)
        this = pd.DataFrame(per_region).mean(axis=1)
        print(f"  {name}: stratified across {len(per_region)} brain regions:",
              ", ".join(sorted(set(region_of))))
    else:
        keep = case | control
        h = g.loc[:, keep].apply(zscore_per_sample)
        this = h.loc[:, case[keep]].mean(axis=1) - h.loc[:, ~case[keep]].mean(axis=1)
    return this.astype(float), g, case, control

deltas, ingested = {}, {}
for name, cfg in PANEL.items():
    d, g, case, control = ingest(name, cfg)
    deltas[name] = d
    ingested[name] = dict(expr=g, case=case, control=control)
    print(name, cfg["acc"],
          f"| case n={int(case.sum())}, control n={int(control.sum())},",
          f"genes={g.shape[0]}, |Δ| median={d.abs().median():.3f}")

assert all(len(d) > 5000 for d in deltas.values()), "disease signatures too small"
json.dump({k: {gene: float(v) for gene, v in d.items()} for k, d in deltas.items()},
          open(RESULTS / "e0_disease_deltas.json", "w"))
