import urllib.request, urllib.error, hashlib
from datetime import datetime, timezone

MANIFEST = []
UA = {"User-Agent": "srcwm-e0-notebook/1.0 (research; issues: github.com/behraj/srcwm)"}

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f, open("/dev/null", "w") as _:
        while True:
            b = f.read(1 << 21)
            if not b: break
            h.update(b)
    return h.hexdigest()

def fetch(url, dest, min_bytes=100, retries=5):
    dest = Path(dest); dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        print(f"  cached: {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        MANIFEST.append(dict(file=str(dest), url=url, sha256=sha256_file(dest),
                             bytes=dest.stat().st_size, cached=True,
                             fetched_utc=datetime.now(timezone.utc).isoformat()))
        return dest
    wait = 2.0
    for attempt in range(1, retries + 1):
        try:
            t0 = time.time()
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=900) as r, open(dest, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk: break
                    f.write(chunk)
            size = dest.stat().st_size
            if size < min_bytes:
                raise IOError(f"download suspiciously small: {size} bytes")
            print(f"  ok: {dest.name} {size/1e6:.1f} MB in {time.time()-t0:.0f}s")
            MANIFEST.append(dict(file=str(dest), url=url, sha256=sha256_file(dest),
                                 bytes=size,
                                 fetched_utc=datetime.now(timezone.utc).isoformat()))
            return dest
        except Exception as e:
            if dest.exists(): dest.unlink()
            if attempt == retries: raise
            time.sleep(wait); wait *= 2

print("fetch utility ready (retries with backoff, SHA-256 manifest)")
