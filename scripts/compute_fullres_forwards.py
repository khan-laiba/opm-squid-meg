#!/usr/bin/env python3
"""Full-resolution cortical lead fields (sample subject, 311,994 valid white-surface vertices)
for G1B, G1C and G2. Resumable: every 20k-source chunk is cached by content hash
(cache/fwd), and the assembled float32 matrices are written to cache/fullres/<job>.npy.

Jobs
----
neuromag4pt_hunold  Neuromag T3 with the 4-point rule (Hunold), BEM 0.33/0.0042/0.33 S/m
opm_hunold          matched OPM array (99 sites), BEM 0.33/0.0042/0.33 S/m
neuromag_bem006     Neuromag T3 (MNE accurate rule), BEM 0.3/0.006/0.3 S/m (MNE default)
opm_bem006          matched OPM array, BEM 0.3/0.006/0.3 S/m
neuromag_bem06      Neuromag T3, BEM 0.3/0.06/0.3 S/m (Goldenholz et al. as printed)

Usage: python scripts/compute_fullres_forwards.py [job ...]   (default: all, in this order)
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import mne  # noqa: E402

from opmsquid import anatomy, forward, neuromag, opm, paths  # noqa: E402

BEMS = {"hunold": (0.33, 0.0042, 0.33), "bem006": (0.3, 0.006, 0.3), "bem06": (0.3, 0.06, 0.3)}
JOBS = {
    "neuromag4pt_hunold": ("T3-4pt", "hunold"),
    "opm_hunold": ("opm", "hunold"),
    "neuromag_bem006": ("T3", "bem006"),
    "opm_bem006": ("opm", "bem006"),
    "neuromag_bem06": ("T3", "bem06"),
}


def array_info(kind: str, subject) -> mne.Info:
    if kind == "opm":
        dig = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
        arr, _ = opm.matched_to_neuromag(neuromag.load_info("T3"), subject.trans, subject.scalp, dig)
        return opm.make_info(arr)
    return neuromag.load_info(kind)


def main(jobs):
    mne.set_log_level("WARNING")
    subject = anatomy.load_sample()
    cortex = anatomy.full_resolution(subject)
    idx = np.flatnonzero(cortex.valid)
    out_dir = paths.CACHE / "fullres"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / "valid_index.npy", idx)
    for job in jobs:
        target = out_dir / f"{job}.npy"
        if target.exists():
            print(f"{job}: exists, skipped")
            continue
        kind, bem_name = JOBS[job]
        t0 = time.time()
        info = array_info(kind, subject)
        gain = forward.chunked_discrete_gain(info, subject.trans, cortex.rr[idx], cortex.nn[idx],
                                             subject.bem_model(BEMS[bem_name]), coil_def=opm.coil_def_file(),
                                             label=job)
        np.save(target, gain)
        with open(out_dir / f"{job}.channels.txt", "w") as fh:
            fh.write("\n".join(info.ch_names))
        print(f"{job}: {gain.shape} in {time.time() - t0:.0f} s -> {target}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or list(JOBS))
