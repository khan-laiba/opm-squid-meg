#!/usr/bin/env python3
"""Fetch the school-aged children of OpenNeuro ds005234 for G3B, the pediatric G4 studies and their MRI quality check.

Two manifests list the files, each with the S3 object version it is downloaded from, its size and its SHA-256:

  configs/school_subjects_manifest.json (45 files, 124,961,312 bytes; snapshot 2.2.0): the children's white, sphere
      and sphere.reg surfaces, aparc annotations and dense MRI scalps; the talairach.xfm and watershed BEM files stored
      in their folders, which the snapshot's file tree shifts by one subject (child A's own are in sub-Z209's folder;
      those in sub-Z213's and sub-Z226's belong to sub-Z214 and sub-Z227: their talairach.xfm serves only
      scripts/study_school_anatomy.py, their BEM files are not used); and children B's and C's own BEMs, from
      sub-Z207's and sub-Z224's folders. Then run scripts/prepare_school_subjects.py.
  configs/school_subjects_qc_manifest.json (6 files, 9,855,337 bytes; the same snapshot): the children's T1 volumes
      and the head masks their scalps were tessellated from (mri/T1.mgz, mri/seghead.mgz), read only by the MRI quality
      check (scripts/study_children_qc.py, which without them reports the children as 'undetermined'). --no-qc leaves
      them out.

Every file goes to data/external/school_subjects/<dest> (OPMSQUID_DATA; src/opmsquid/paths.py). It is downloaded from
its exact S3 object version (https://s3.amazonaws.com/<bucket>/<key>?versionId=<version>; the quality-check manifest
gives the object as s3://<bucket>/<key> with its version_id) into <dest>.part, which is checked against the manifest's
size and SHA-256 and only then moved into place: a download that differs, or breaks off, is deleted and the script
stops. Files already present and intact are kept, so a second run downloads nothing.

  .venv/bin/python scripts/fetch_school_subjects.py            # download what is missing (both manifests)
  .venv/bin/python scripts/fetch_school_subjects.py --no-qc    # the same without the quality-check files
  .venv/bin/python scripts/fetch_school_subjects.py --check    # check only, download nothing; exit status 1 if a
                                                               # file is missing or changed (--no-qc: surfaces only)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from opmsquid import anatomy, paths  # noqa: E402

MANIFEST = ROOT / "configs" / "school_subjects_manifest.json"
QC_MANIFEST = ROOT / "configs" / "school_subjects_qc_manifest.json"
MANIFESTS = (("surfaces", MANIFEST), ("qc", QC_MANIFEST))  # checked and fetched in this order; --no-qc drops "qc"
LABELS = {"surfaces": "surfaces, BEMs and transforms", "qc": "MRI quality check: T1 volumes and head masks"}
S3_HTTPS = "https://s3.amazonaws.com/"  # path-style address of s3://<bucket>/<key> (the bucket's name has dots)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def intact(path: Path, entry: dict) -> bool:
    return path.is_file() and path.stat().st_size == entry["size"] and sha256(path) == entry["sha256"]


def entry(raw: dict, where: str = "manifest") -> dict:
    """One listed file as fetched here: dest (relative to the subjects folder), size, sha256, version, source and url,
    the address of that exact S3 object version. The surfaces' manifest gives the object's https address and `version`,
    the quality check's its s3:// URI and `version_id`. ValueError for an entry that cannot be fetched and checked
    exactly (no object version, a destination outside the folder, a malformed size or checksum, an address that is not
    https or s3://)."""
    lacking = {"dest", "size", "sha256", "source"} - set(raw)
    if lacking:
        raise ValueError(f"{where}: no {', '.join(sorted(lacking))}")
    dest = PurePosixPath(raw["dest"]) if isinstance(raw["dest"], str) else None
    if dest is None or dest.is_absolute() or not dest.parts or ".." in dest.parts:
        raise ValueError(f"{where}: dest {raw['dest']!r} is not a relative path inside the subjects folder")
    size = raw["size"]
    if not isinstance(size, int) or isinstance(size, bool) or size < 0:
        raise ValueError(f"{where}: size {size!r} is not a byte count")
    if not isinstance(raw["sha256"], str) or not re.fullmatch(r"[0-9a-f]{64}", raw["sha256"]):
        raise ValueError(f"{where}: sha256 {raw['sha256']!r} is not 64 lower-case hexadecimal digits")
    version = raw.get("version", raw.get("version_id"))
    if not isinstance(version, str) or not version:
        raise ValueError(f"{where}: no S3 object version (version or version_id)")
    source = raw["source"]
    if isinstance(source, str) and source.startswith("s3://"):
        source = S3_HTTPS + source[len("s3://"):]
    parts = urllib.parse.urlsplit(source) if isinstance(source, str) else None
    if parts is None or parts.scheme != "https" or not parts.netloc or parts.query or parts.fragment:
        raise ValueError(f"{where}: source {raw['source']!r} is not an https or s3:// address without a query")
    return {"dest": str(dest), "size": size, "sha256": raw["sha256"], "version": version, "source": source,
            "url": f"{source}?versionId={urllib.parse.quote(version, safe='')}"}


def load(path: Path) -> dict:
    """The manifest at `path` with its files as entry() gives them. A file count (n_files) or byte total (total_bytes)
    that the manifest states must be that of its files, and no destination may be listed twice (ValueError)."""
    path = Path(path)
    manifest = json.loads(path.read_text())
    if not isinstance(manifest.get("files"), list):
        raise ValueError(f"{path.name}: no list of files")
    files = [entry(e, f"{path.name}, file {i + 1}") for i, e in enumerate(manifest["files"])]
    for key, value in (("n_files", len(files)), ("total_bytes", sum(e["size"] for e in files))):
        if key in manifest and manifest[key] != value:
            raise ValueError(f"{path.name}: {key} is {manifest[key]}, its files give {value}")
    if len({e["dest"] for e in files}) != len(files):
        raise ValueError(f"{path.name}: a destination is listed twice")
    return {**manifest, "files": files}


def fetch(entry: dict, dest: Path) -> None:
    """Download `entry` (its exact S3 object version) to `dest` through dest.part, moved into place only when its size
    and SHA-256 are the manifest's; when they are not, or the download fails, the .part file is deleted."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    try:
        with urllib.request.urlopen(entry["url"], timeout=120) as r, open(part, "wb") as fh:
            while block := r.read(1 << 20):
                fh.write(block)
        if not intact(part, entry):
            raise RuntimeError(f"{entry['dest']}: size or SHA-256 differs from the manifest after download")
    except BaseException:
        part.unlink(missing_ok=True)
        raise
    part.replace(dest)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Fetch, or check, the files of the school-aged children listed in configs/"
                                             "school_subjects_manifest.json and configs/school_subjects_qc_manifest.json.")
    ap.add_argument("--check", action="store_true",
                    help="check only, download nothing; exit status 1 if a file is missing or changed")
    ap.add_argument("--no-qc", action="store_true",
                    help="leave out the MRI quality-check files (configs/school_subjects_qc_manifest.json)")
    args = ap.parse_args(argv)
    manifests = [(key, Path(path), load(path)) for key, path in MANIFESTS if not (args.no_qc and key == "qc")]
    dests = [e["dest"] for _, _, m in manifests for e in m["files"]]
    if len(set(dests)) != len(dests):
        raise ValueError("a destination is listed in more than one manifest")
    root = paths.EXTERNAL / anatomy.SCHOOL_SUBJECTS
    first = manifests[0][2]
    print(" ".join(str(first[k]) for k in ("dataset", "snapshot") if k in first) + f" in {root}:")
    todo = []
    for key, path, m in manifests:
        missing = [e for e in m["files"] if not intact(root / e["dest"], e)]
        n = len(m["files"])
        print(f"  {path.name} ({LABELS.get(key, key)}): {n - len(missing)} of {n} files present and intact")
        todo += missing
    if args.check or not todo:
        for e in todo:
            print(f"  missing or changed: {e['dest']} ({e['size']:,d} bytes)")
        return 1 if todo else 0
    print(f"downloading {len(todo)} files, {sum(e['size'] for e in todo):,d} bytes "
          f"(cite: {first.get('citation', first.get('dataset'))}; licence: {first.get('license')})", flush=True)
    for e in todo:
        fetch(e, root / e["dest"])
        print(f"  {e['dest']}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
