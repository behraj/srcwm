print("LINCS L1000 metadata (GSE92742, Phase 1) ...")
for f in LINCS_FILES:
    fetch(LINCS_BASE + f, DATA / "lincs" / f)

checksums = {}
cl = (DATA / "lincs" / "GSE92742_SHA512SUMS.txt.gz").read_bytes()
for line in gzip.decompress(cl).decode().splitlines():
    parts = line.split()
    if len(parts) >= 2:
        checksums[parts[1]] = parts[0]
print("published SHA512SUMS entries:", len(checksums))

n_verified = 0
for m in MANIFEST:
    base = m["file"].split("/")[-1]
    if base in checksums:
        h = hashlib.sha512(open(m["file"], "rb").read()).hexdigest()
        assert h == checksums[base], f"LINCS SHA512 mismatch: {base}"
        n_verified += 1
        print(f"  LINCS SHA512 verified: {base} "
              f"{'(cache hit)' if m.get('cached') else '(fresh download)'}")
assert n_verified >= 1, "no LINCS file matched a published SHA512 entry — verification did not run"
print(f"LINCS SHA512 verification: {n_verified} file(s) match published sums")

print("GEO calibration series + platform annotations ...")
for name, cfg in PANEL.items():
    acc, gpl = cfg["acc"], cfg["gpl"]
    fetch(f"{geo_dir(acc)}/matrix/{acc}_series_matrix.txt.gz",
          DATA / "geo" / f"{acc}_series_matrix.txt.gz")
    fetch(f"{gpl_dir(gpl)}/annot/{gpl}.annot.gz", DATA / "geo" / f"{gpl}.annot.gz")

print("STRING v12.0 (9606) ...")
for url, name in STRING_FILES:
    fetch(url, DATA / "string" / name)
print("all static downloads complete")
