import urllib.request

def chembl_get(path_and_query):
    url = CHEMBL_API + path_and_query
    req = urllib.request.Request(url, headers=UA)
    for i in range(4):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=180).read().decode())
        except Exception:
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"ChEMBL request failed: {url}")

ind_path = DATA / "snapshots" / "chembl_drug_indication_snapshot.json"
mol_path = DATA / "snapshots" / "chembl_molecules.json"
if ind_path.exists() and mol_path.exists():
    ind = pd.read_json(ind_path, lines=True)
    mol = pd.read_json(mol_path, lines=True)
    print("ChEMBL snapshots loaded from data/snapshots (cached):", ind.shape, mol.shape)
else:
    rows, offset = [], 0
    while True:
        d = chembl_get(f"/drug_indication.json?limit={CHEMBL_PAGE}&offset={offset}")
        page = d["drug_indications"]
        rows.extend(page)
        if len(page) < CHEMBL_PAGE:
            break
        offset += CHEMBL_PAGE
        if offset > 200_000:
            raise RuntimeError("drug_indication pagination runaway")
    ind = pd.DataFrame(rows)
    ind["efo_term_l"] = ind.efo_term.fillna("").str.lower()
    ind["max_phase_for_ind_f"] = pd.to_numeric(ind.max_phase_for_ind, errors="coerce")
    print("drug_indication rows:", len(ind), "| unique drugs:", ind.molecule_chembl_id.nunique())

    want_ids = sorted(ind.molecule_chembl_id.unique())
    print("fetching pref_names for", len(want_ids), "indicated drugs (chunked __in filter)")
    names = []
    for s in range(0, len(want_ids), 200):
        chunk = want_ids[s:s + 200]
        d = chembl_get("/molecule.json?molecule_chembl_id__in=" + ",".join(chunk))
        for m in d["molecules"]:
            names.append(dict(molecule_chembl_id=m["molecule_chembl_id"],
                              pref_name=(m.get("pref_name") or "").lower(),
                              max_phase=m.get("max_phase")))
        time.sleep(0.2)
    mol = pd.DataFrame(names).drop_duplicates("molecule_chembl_id")

    ind_path.parent.mkdir(parents=True, exist_ok=True)
    ind.to_json(ind_path, orient="records", lines=True, force_ascii=False)
    mol.to_json(mol_path, orient="records", lines=True)
    print("ChEMBL snapshot frozen:", ind.shape, mol.shape)
