import os, io, json, gzip, time, re, hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timezone

SMOKE = True

DATA = Path("data"); RESULTS = Path("results"); FIGS = Path("figures")
for p in (DATA, DATA/"geo", DATA/"lincs", DATA/"string", DATA/"snapshots", RESULTS, FIGS):
    p.mkdir(parents=True, exist_ok=True)

SEED = 0
rng = np.random.default_rng(SEED)

CENSUS_VERSION = "2025-11-08"

LINCS_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/"
LINCS_FILES = [
    "GSE92742_Broad_LINCS_gene_info.txt.gz",
    "GSE92742_Broad_LINCS_pert_info.txt.gz",
    "GSE92742_Broad_LINCS_inst_info.txt.gz",
    "GSE92742_Broad_LINCS_sig_info.txt.gz",
    "GSE92742_Broad_LINCS_cell_info.txt.gz",
    "GSE92742_SHA512SUMS.txt.gz",
]
LINCS_LEVEL5 = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"

def _geo_record_dir(prefix, n):
    grp = "nnn" if n < 1000 else f"{n // 1000}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/{prefix}{grp}"

def geo_dir(acc):
    return f"{_geo_record_dir('series/GSE', int(acc[3:]))}/{acc}"
def gpl_dir(gpl):
    return f"{_geo_record_dir('platforms/GPL', int(gpl[3:]))}/{gpl}"

PANEL = {
    "AD": dict(
        acc="GSE5281", gpl="GPL570", tissue="brain (laser-captured neurons, 7 regions)",
        case_terms=["alzheimer"],
        control_terms=["disease state: normal"],
        region_line_contains="organ region",
        stratify="region",
        census_tissue_general="brain",
        expected_cell_types=["microglial cell", "astrocyte", "neuron",
                             "glutamatergic neuron"],
    ),
    "RA": dict(
        acc="GSE55235", gpl="GPL96", tissue="synovium",
        case_terms=["disease state: rheumatoid arthritis"],
        control_terms=["disease state: healthy control"],
        stratify=None,
        census_tissue_general="blood",
        expected_cell_types=["CD14-positive monocyte", "classical monocyte",
                             "T cell", "B cell"],
        note=("RA bulk series is synovial tissue; a synovial single-cell reference "
              "is not available from a keyless source (the Census has none), so the "
              "lifting demo in this notebook uses blood as the reference tissue and "
              "reports this limitation explicitly. Full E3 needs a synovial atlas."),
    ),
    "IPF": dict(
        acc="GSE24206", gpl="GPL570", tissue="lung",
        case_terms=["idiopathic pulmonary fibrosis"],
        control_terms=["phenotype: healthy"],
        stratify=None,
        census_tissue_general="lung",
        expected_cell_types=["pulmonary alveolar type 2 cell", "macrophage",
                             "alveolar macrophage"],
    ),
}

STRING_RELEASE = "v12.0"
STRING_BASE = "https://stringdb-downloads.org/download"
STRING_FILES = [
    (f"{STRING_BASE}/protein.links.{STRING_RELEASE}/9606.protein.links.{STRING_RELEASE}.txt.gz",
     "9606.protein.links.v12.0.txt.gz"),
    (f"{STRING_BASE}/protein.info.{STRING_RELEASE}/9606.protein.info.{STRING_RELEASE}.txt.gz",
     "9606.protein.info.v12.0.txt.gz"),
]

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_PAGE = 1000
INDICATION_TERMS = {
    "AD":  ["alzheimer"],
    "RA":  ["rheumatoid"],
    "IPF": ["pulmonary fibrosis"],
}

PRIMARY_METRICS = ["P@25", "NDB", "dose-consistency",
                   "program-rescue coverage", "safety-margin diversity"]
E5_DECISION_RULE = ("SRCWM beats M-LIN on >= 3 of 5 LOGO folds at p<0.05 (Holm); "
                    "otherwise the conclusion is negative and the emphasis shifts to "
                    "the E6 regime analysis. Binding; no post-hoc switch.")

TOP_K = 500          # official up/down sets for signature construction
PCA_PROGRAMS = 50    # latent program dimension
ROLLOUT_T = [1, 2, 3]
BOOTSTRAP_B = 2000   # SMOKE; pre-registered value is 10000 for E5

CONFIG = dict(SMOKE=SMOKE, seed=SEED, census_version=CENSUS_VERSION,
              panel={k: {kk: (sorted(vv) if isinstance(vv, (list, set)) else vv)
                         for kk, vv in v.items()} for k, v in PANEL.items()},
              string_release=STRING_RELEASE,
              primary_metrics=PRIMARY_METRICS,
              e5_decision_rule=E5_DECISION_RULE,
              top_k=TOP_K, pca_programs=PCA_PROGRAMS,
              rollout_T=ROLLOUT_T, bootstrap_B=BOOTSTRAP_B,
              created_utc=datetime.now(timezone.utc).isoformat())
CONFIG_HASH = hashlib.sha256(json.dumps(CONFIG, sort_keys=True, default=str)
                             .encode()).hexdigest()[:16]
print("config hash (E10 experiment key prefix):", CONFIG_HASH)
print("SMOKE:", SMOKE)
