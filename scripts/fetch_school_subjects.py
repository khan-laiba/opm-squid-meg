#!/usr/bin/env python3
"""Fetch the school-aged children of OpenNeuro ds005234 for G3B and the pediatric G4 studies.

The 45 files of `configs/school_subjects_manifest.json` (124,961,312 bytes; snapshot 2.2.0) are
downloaded from the S3 object versions listed there into data/external/school_subjects/<dest> and
checked against their size and SHA-256; files already present and intact are kept. They are the
children's white, sphere and sphere.reg surfaces, aparc annotations and dense MRI scalps; the
talairach.xfm and watershed BEM files stored in their folders, which the snapshot's file tree
shifts by one subject (child A's own are in sub-Z209's folder; those in sub-Z213's and sub-Z226's
belong to sub-Z214 and sub-Z227: their talairach.xfm serves only scripts/study_school_anatomy.py,
their BEM files are not used); and children B's and C's own BEMs, from sub-Z207's and sub-Z224's
folders. Then run scripts/prepare_school_subjects.py.

  .venv/bin/python scripts/fetch_school_subjects.py           # download what is missing
  .venv/bin/python scripts/fetch_school_subjects.py --check   # check only, download nothing
"""
from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from opmsquid import anatomy, paths  # noqa: E402

MANIFEST = ROOT / "configs" / "school_subjects_manifest.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def intact(path: Path, entry: dict) -> bool:
    return path.is_file() and path.stat().st_size == entry["size"] and sha256(path) == entry["sha256"]


def url(entry: dict) -> str:
    return f"{entry['source']}?versionId={entry['version']}"


def fetch(entry: dict, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    with urllib.request.urlopen(url(entry), timeout=120) as r, open(part, "wb") as fh:
        while block := r.read(1 << 20):
            fh.write(block)
    if not intact(part, entry):
        part.unlink()
        raise RuntimeError(f"{entry['dest']}: size or SHA-256 differs from the manifest after download")
    part.replace(dest)


def main(argv: list[str]) -> int:
    check_only = "--check" in argv
    manifest = json.loads(MANIFEST.read_text())
    root = paths.EXTERNAL / anatomy.SCHOOL_SUBJECTS
    missing = [e for e in manifest["files"] if not intact(root / e["dest"], e)]
    print(f"{manifest['dataset']} {manifest['snapshot']}: {manifest['n_files'] - len(missing)} of {manifest['n_files']} "
          f"files present and intact in {root}")
    if check_only or not missing:
        for e in missing:
            print(f"  missing or changed: {e['dest']} ({e['size']:,d} bytes)")
        return 1 if missing else 0
    print(f"downloading {len(missing)} files, {sum(e['size'] for e in missing):,d} bytes "
          f"(cite: {manifest['citation']}; licence: {manifest['license']})", flush=True)
    for e in missing:
        fetch(e, root / e["dest"])
        print(f"  {e['dest']}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
