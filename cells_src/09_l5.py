gctx_path = DATA / "lincs" / LINCS_LEVEL5
if SMOKE:
    print("SMOKE profile: the 21 GB LINCS Level-5 GCTX is NOT downloaded.")
    print("The streaming reader is validated below on a synthetic GCTX with the")
    print("identical HDF5 layout; set SMOKE=False to stream the real file.")
else:
    print("FULL: downloading", LINCS_LEVEL5, "(~21 GB) ...")
    fetch(LINCS_BASE + LINCS_LEVEL5, gctx_path)
    expected = checksums.get(LINCS_LEVEL5)
    if expected:
        h = hashlib.sha512(gctx_path.read_bytes()).hexdigest()
        assert h == expected, "LINCS Level-5 SHA512 mismatch"
        print("LINCS Level-5 SHA512 verified.")
