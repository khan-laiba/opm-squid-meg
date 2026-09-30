"""Result writing helpers: strict JSON (NaN/inf -> null, numpy -> Python) and provenance stamps."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np

from . import paths


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


CODE_PATHS = ("src", "scripts", "configs", "tests", "requirements.txt")


def git_commit() -> str:
    """Current commit, with '+dirty' if code or configuration (``CODE_PATHS``, untracked files
    included) differ from it. Generated outputs (results/, docs) do not count: a script
    regenerating its own tracked figures must not mark its provenance dirty."""
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=paths.ROOT, capture_output=True,
                             text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", *CODE_PATHS], cwd=paths.ROOT,
                               capture_output=True, text=True, check=True).stdout.strip()
        return sha + ("+dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def write_json(obj: dict, path: Path) -> None:
    import mne

    obj = dict(obj)
    obj.setdefault("provenance", dict(commit=git_commit(), mne_version=mne.__version__, numpy_version=np.__version__))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(json_safe(obj), fh, indent=2, allow_nan=False)
