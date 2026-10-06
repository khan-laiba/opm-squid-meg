"""Locations of external data, caches and results.

External data (MNE sample data, anatomy) and caches are never committed. Their locations can
be overridden with the environment variables ``OPMSQUID_DATA`` and ``OPMSQUID_CACHE``.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("OPMSQUID_DATA", ROOT / "data"))
EXTERNAL = DATA / "external"
SAMPLE_DATA = EXTERNAL / "MNE-sample-data"  # mne.datasets.sample, MNE 1.13.2
SAMPLE_MEG = SAMPLE_DATA / "MEG" / "sample"
SUBJECTS_DIR = SAMPLE_DATA / "subjects"
CACHE = Path(os.environ.get("OPMSQUID_CACHE", ROOT / "cache"))
RESULTS = ROOT / "results"
CONFIGS = ROOT / "configs"
LEGACY = ROOT / "legacy"


def require(path: Path, what: str) -> Path:
    """Return `path` if it exists; otherwise fail with a message naming the missing input."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{what} not found at {path} (see README.md, 'External data')")
    return path
