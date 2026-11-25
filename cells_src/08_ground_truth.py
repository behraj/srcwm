name_by_chembl = dict(zip(mol.molecule_chembl_id, mol.pref_name))
ground_truth = {}
rows = []
for disease, terms in INDICATION_TERMS.items():
    m = ind[ind.efo_term_l.apply(lambda s: any(t in s for t in terms))]
    ph2 = m[m.max_phase_for_ind_f >= 2.0]
    chembl_ids = sorted(set(ph2.molecule_chembl_id))
    names = sorted({n for n in ph2.molecule_chembl_id.map(name_by_chembl)
                    if isinstance(n, str) and n})
    lincs_names = lincs_drug_names
    in_l1000 = sorted({n for n in names if n in lincs_names})
    ground_truth[disease] = dict(chembl_ids=chembl_ids, drug_names=names,
                                 in_l1000=in_l1000)
    rows.append((disease, m.shape[0], ph2.shape[0], len(names), len(in_l1000)))
gt = pd.DataFrame(rows, columns=["disease", "indications (all phases)",
                                 "indications phase>=2", "named drugs (phase>=2)",
                                 "present in L1000"])
gt
