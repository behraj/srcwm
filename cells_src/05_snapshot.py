from datetime import datetime, timezone

def gpl_dir_for_md(gpl):
    return gpl_dir(gpl)

lines = []
lines.append("# SRCWM frozen data snapshot\n")
lines.append(f"- Created: {datetime.now(timezone.utc).isoformat()}")
lines.append(f"- Notebook: `srcwm_e0_calibration.ipynb`")
lines.append(f"- E10 experiment key (config hash): `{CONFIG_HASH}`")
lines.append(f"- Profile: **{'SMOKE' if SMOKE else 'FULL'}**\n")
lines.append("Every component below was fetched through this notebook directly from the")
lines.append("primary source. Local SHA-256 checksums are recomputed after download; LINCS")
lines.append("files are additionally verified against the published GSE92742 SHA512SUMS.\n")
lines.append("## Components\n")
lines.append("| Component | Source | Notes |")
lines.append("|---|---|---|")
lines.append(f"| LINCS L1000 metadata | {LINCS_BASE} | GSE92742 Phase 1 |")
if not SMOKE:
    lines.append(f"| LINCS L1000 Level 5 signatures | {LINCS_BASE + LINCS_LEVEL5} | 21 GB, "
                 "positions streamed, never fully materialised |")
lines.append("| STRING PPI protein.links.v12.0 (9606) | STRING v12.0 | frozen |")
lines.append("| GEO bulk series | " + ", ".join(sorted({c['acc'] for c in PANEL.values()})) +
             " | series matrices + NCBI annot files |")
lines.append(f"| ChEMBL indications | {CHEMBL_API}/drug_indication + molecule | paged "
             "API snapshot, frozen to data/snapshots/ |")
lines.append(f"| Single-cell reference | CZ CELLxGENE Census `{CENSUS_VERSION}` (LTS) | "
             "queried remotely; slice hashes in results |\n")
lines.append("## Local files\n")
lines.append("| file | bytes | sha256 (first 16) |")
lines.append("|---|---|---|")
for m in MANIFEST:
    lines.append(f"| {m['file']} | {m['bytes']:,} | {m['sha256'][:16]}… |")
(DATA / "SNAPSHOT.md").write_text("\n".join(lines) + "\n")
print((DATA / "SNAPSHOT.md").read_text()[:1400])
