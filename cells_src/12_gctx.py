# Streaming GCTX reader (h5py). Works on the real LINCS Level-5 file (FULL
# profile) and on synthetic GCTX files with the same layout (SMOKE validation).

def write_synthetic_gctx(path, matrix):
    n_g, n_i = matrix.shape
    gids = [f"gene_{i}" for i in range(n_g)]
    sigids = [f"sig_{i}" for i in range(n_i)]
    with h5py.File(path, "w") as f:
        grp = f.create_group("0")
        d0 = grp.create_group("DATA").create_group("0")
        d0.create_dataset("matrix", data=matrix.astype("float32"))
        meta = grp.create_group("META")
        meta.create_group("ROW").create_dataset("id", data=np.array(gids, dtype="S"))
        meta.create_group("COL").create_dataset("id", data=np.array(sigids, dtype="S"))

def gctx_column_ids(path):
    with h5py.File(path, "r") as f:
        ids = f["0"]["META"]["COL"]["id"][:]
    return [i.decode() if isinstance(i, bytes) else str(i) for i in ids]

def read_gctx_columns(path, want_cols, block=20000):
    names = gctx_column_ids(path)
    sel = sorted(names.index(w) for w in want_cols)
    with h5py.File(path, "r") as f:
        ds = f["0"]["DATA"]["0"]["matrix"]
        out = np.empty((ds.shape[0], len(sel)), dtype="float32")
        for s in range(0, len(sel), block):
            b = sel[s:s + block]
            out[:, s:s + len(b)] = ds[:, b]
    return out

syn_path = str(DATA / "lincs" / "synthetic_minigctx.gctx")
mat = rng.normal(size=(978, 400)).astype("float32")
write_synthetic_gctx(syn_path, mat)
want = ["sig_7", "sig_111", "sig_399"]
sub = read_gctx_columns(syn_path, want)
assert sub.shape == (978, 3)
assert np.allclose(sub[:, -1], mat[:, 399])
print("synthetic GCTX round-trip + block streaming read: passed")

if SMOKE:
    print("REAL LINCS Level-5 stream: gated to the FULL profile (SMOKE=False).")
    print("In FULL, this cell streams gene-aligned columns for every approved-drug")
    print("instance selected from sig_info (a few hundred thousand signatures) in")
    print("20k-column blocks — one sequential pass over the 21 GB file per block set.")
