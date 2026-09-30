#!/usr/bin/env python3
"""Full-resolution cortical lead fields (sample subject, 311,994 valid white-surface vertices)
for G1B, G1C, G2 and G4 (job table in opmsquid.fullres). Resumable: every 20k-source chunk is
cached by content hash (cache/fwd); the assembled float32 matrices go to cache/fullres/<job>.npy
with a fingerprint sidecar (<job>.meta.json). A job whose stored fingerprint matches the current
array and BEM is kept; a stale or unfingerprinted one is recomputed.

Jobs
----
neuromag4pt_hunold  Neuromag T3 with the 4-point rule (Hunold), BEM 0.33/0.0042/0.33 S/m
opm_hunold          matched OPM array, BEM 0.33/0.0042/0.33 S/m
neuromag_bem006     Neuromag T3 (MNE accurate rule), BEM 0.3/0.006/0.3 S/m (MNE default)
opm_bem006          matched OPM array, BEM 0.3/0.006/0.3 S/m
neuromag_bem06      Neuromag T3, BEM 0.3/0.06/0.3 S/m (Goldenholz et al. as printed)
opm204_bem006       dense OPM array, 204 sites (G2 channel-budget control), BEM 0.3/0.006/0.3 S/m
opm_dense_bem006    dense OPM array (G2 full system), BEM 0.3/0.006/0.3 S/m

Usage: python scripts/compute_fullres_forwards.py [--force] [job ...]   (default: all, in this order)
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import mne  # noqa: E402

from opmsquid import anatomy, fullres  # noqa: E402

JOBS = fullres.JOBS


def main(jobs, force=False):
    mne.set_log_level("WARNING")
    subject = anatomy.load_sample()
    cortex = anatomy.full_resolution(subject)
    for job in jobs:
        fullres.compute(job, subject, cortex, force=force, log=lambda m: print(m, flush=True))


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--force"]
    main(args or list(JOBS), force="--force" in sys.argv)
