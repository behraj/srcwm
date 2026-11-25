def read_gz_csv(path):
    with gzip.open(path, "rt") as f:
        return pd.read_csv(f, sep="\t", dtype=str)

gene_info = read_gz_csv(DATA / "lincs" / "GSE92742_Broad_LINCS_gene_info.txt.gz")
gene_info["pr_gene_id"] = gene_info["pr_gene_id"].astype(int)
gene_info["pr_gene_symbol_u"] = gene_info.pr_gene_symbol.str.upper()
pert_info = read_gz_csv(DATA / "lincs" / "GSE92742_Broad_LINCS_pert_info.txt.gz")
inst_info = read_gz_csv(DATA / "lincs" / "GSE92742_Broad_LINCS_inst_info.txt.gz")
sig_info = read_gz_csv(DATA / "lincs" / "GSE92742_Broad_LINCS_sig_info.txt.gz")
cell_info = read_gz_csv(DATA / "lincs" / "GSE92742_Broad_LINCS_cell_info.txt.gz")

print("gene_info:", gene_info.shape, "| landmark genes:", (gene_info.pr_is_lm == "1").sum())
print("pert_info:", pert_info.shape, "| pert types:", pert_info.pert_type.value_counts().to_dict())
print("inst_info:", inst_info.shape)
print("sig_info:", sig_info.shape, "| compound signature rows:", (sig_info.pert_type == "trt_cp").sum())
pd.set_option("display.width", 160)
trt_cp = sig_info[sig_info.pert_type == "trt_cp"]
print("top cell lines by compound signature count:")
print(trt_cp.groupby("cell_id").size().sort_values(ascending=False).head(9).to_string())
print("distinct recombinant drugs (trt_cp):", trt_cp.pert_iname.str.lower().nunique())

N_TRT_CP_SIGS = int(len(trt_cp))
lincs_drug_names = set(pert_info.pert_iname.str.lower())
del trt_cp, inst_info, sig_info, cell_info
import gc; gc.collect()
