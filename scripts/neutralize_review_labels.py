#!/usr/bin/env python3
"""Replace the review-process labels in four stored revision results with the neutral descriptions their scripts now write.

The revision's checks first recorded, in their output metadata, which review comments they answered (verbatim quotations
and the reviewers' labels). Their scripts now write neutral descriptions of the analyses instead (a `purpose` field and a
neutral `status`). This script brings the stored files in line without re-running anything: it touches only those
metadata fields, keeps every number and the recorded provenance (the commit that computed the results) unchanged, and
records the commit that revised the metadata in `metadata_revised_at_commit`. It checks that nothing else changed.

Files: results/g2_covariance_validation/covariance_validation.json, results/g2_noise_sensitivity/
noise_sensitivity_summary.json, results/g3b_children_qc/children_qc.json, results/g2/g2_depth_bins.json.
Usage: .venv/bin/python scripts/neutralize_review_labels.py
"""
from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from opmsquid import io  # noqa: E402

OLD_NOTE = "thresholds 8 and 10 mm from the referee's request "
NEW_NOTE = "thresholds 8 and 10 mm as in the shallow-cortex question "
OLD_CONTRAST = "is the referee's contrast."
NEW_CONTRAST = "is the contrast it is compared with."


def _module(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _replace_key(d: dict, old: str, new: str, value) -> dict:
    """``d`` with key ``old`` replaced by ``new`` (same position) holding ``value``."""
    if old not in d:
        raise KeyError(old)
    return {(new if k == old else k): (value if k == old else v) for k, v in d.items()}


def main():
    cov = _module("study_covariance_validation")
    ns = _module("study_noise_sensitivity")
    qc = _module("study_children_qc")
    db = _module("study_g2_depth_bins")
    edits = {
        "results/g2_covariance_validation/covariance_validation.json":
            lambda d: _replace_key(d, "requests", "purpose", cov.PURPOSE),
        "results/g2_noise_sensitivity/noise_sensitivity_summary.json":
            lambda d: _contrast(_replace_key(dict(d, status=ns.STATUS), "referee_requests", "purpose", ns.PURPOSE)),
        "results/g3b_children_qc/children_qc.json":
            lambda d: _note(_replace_key(dict(d, status=qc.STATUS), "requests", "purpose", qc.PURPOSE)),
        "results/g2/g2_depth_bins.json": lambda d: dict(d, status=db.STATUS),
    }
    for rel, fn in edits.items():
        path = ROOT / rel
        before = json.loads(path.read_text())
        after = fn(copy.deepcopy(before))
        after["metadata_revised_at_commit"] = io.RUN_COMMIT
        # nothing but the named metadata may differ
        changed = ({"status", "purpose", "requests", "referee_requests", "metadata_revised_at_commit"}
                   | ({"parameters"} if "parameters" in before else set())
                   | ({"coloured"} if "coloured" in before else set()))
        for k in set(before) | set(after):
            if k not in changed and before.get(k) != after.get(k):
                raise ValueError(f"{rel}: {k} would change")
        if "coloured" in before:
            b = {k: v for k, v in before["coloured"].items() if k != "method"}
            a = {k: v for k, v in after["coloured"].items() if k != "method"}
            if a != b:
                raise ValueError(f"{rel}: the coloured results other than their method text would change")
        if "parameters" in before:
            b = {k: v for k, v in before["parameters"].items() if k != "note"}
            a = {k: v for k, v in after["parameters"].items() if k != "note"}
            if a != b:
                raise ValueError(f"{rel}: parameters other than the note would change")
        text = json.dumps(after)
        for label in ("Fable", "Codex", "fable", "codex", "referee"):
            if label in text:
                raise ValueError(f"{rel}: '{label}' is still present")
        io.write_json(after, path)  # the results' own format; the stored provenance is kept
        print(f"{rel}: metadata revised")


def _contrast(d: dict) -> dict:
    c = dict(d["coloured"])
    if OLD_CONTRAST not in c.get("method", ""):
        raise ValueError("noise_sensitivity_summary.json: the coloured method text has changed; update OLD_CONTRAST")
    c["method"] = c["method"].replace(OLD_CONTRAST, NEW_CONTRAST)
    return dict(d, coloured=c)


def _note(d: dict) -> dict:
    p = dict(d["parameters"])
    if OLD_NOTE not in p.get("note", ""):
        raise ValueError("children_qc.json: the parameters' note has changed; update OLD_NOTE")
    p["note"] = p["note"].replace(OLD_NOTE, NEW_NOTE)
    return dict(d, parameters=p)


if __name__ == "__main__":
    main()
