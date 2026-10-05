#!/usr/bin/env python3
"""Full-precision known-topography detectability of the adult's targets (companion of results/g2/g2_targets.csv).

results/g2/g2_targets.csv prints detectability to four decimals, so a few target-level orderings cannot be read from it
(e.g. a target printed with the same value for an OPM array and Neuromag). This script rebuilds the adult comparison's
noise model in its primary band exactly as scripts/g2_band_sensitivity.py does for 1-40 Hz (same study, arrays, gains,
brain-noise calibration and room field; nothing new is assumed) and writes, for every target, the detectability of
Neuromag (306 channels), the site-matched and the dense OPM array under each noise condition at full precision, the
paired ratios OPM/Neuromag and explicit "OPM ahead" flags (ratio > 1 at full precision).

Checks before writing: every value rounded to four decimals equals the stored table's; the median log2 ratio and the
share of targets with the OPM ahead equal those of results/g2/g2_summary.json (primary.oracle) for every array and
condition.

Output: results/g2/g2_targets_metrics.csv
Usage: .venv/bin/python scripts/export_g2_target_metrics.py
"""
from __future__ import annotations

import csv
import json
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import mne  # noqa: E402

import g2_adult_comparison as G  # noqa: E402
from opmsquid import background, g2, io, noise  # noqa: E402

OUT = ROOT / "results" / "g2" / "g2_targets_metrics.csv"
TABLE = ROOT / "results" / "g2" / "g2_targets.csv"
SUMMARY = ROOT / "results" / "g2" / "g2_summary.json"
ARRAYS = {"squid": "neuromag", "opm_matched": "site_matched", "opm_dense": "dense"}


def main():
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    conds = cfg["conditions"]["all"]
    st = G.Study(cfg)
    bads = cfg["sensors"]["bads"]
    arrays = g2.build_arrays(st.subject, st.dig)
    squid = arrays["squid"]
    good = ~np.isin(squid.info.ch_names, bads)
    grads = good & (squid.kinds == "grad")
    b = cfg["band"]
    filt = noise.AnalysisFilter(fs=st.dig["sfreq"], l_freq=b["l_freq_hz"], h_freq=b["h_freq_hz"], order=b["order"])
    meas = g2.measured_noise(squid.info, filt, bads)
    env = meas["environment"]
    res = {}
    for name in ARRAYS:
        a = arrays[name]
        g_t, g_g = st.gains(a)
        if name == "squid":
            bs = background.calibrate(st.unit_brain(g_g), grads, float(np.nanmedian(meas["brain"][grads])))
        nz = g2.array_noise(a, g_g, st.src.grid_area, bs, env, filt.enbw(), st.opm_asd)
        res[name] = G.evaluate_all(g_t * st.q, a, nz, conds)

    # the four decimals of the stored table
    rows = io.read_csv(TABLE)
    if len(rows) != st.nt:
        raise ValueError(f"{TABLE.name}: {len(rows)} rows for {st.nt} targets")
    for i, v in enumerate(st.src.target):
        if (int(rows[i]["hemi"]), int(rows[i]["vertno"])) != (int(st.cortex.hemi[v]), int(st.cortex.vertno[v])):
            raise ValueError(f"{TABLE.name}: row {i} is not target {i}")
    worst = 0
    for name in ARRAYS:
        cs = "combined" if name == "squid" else "opm"
        for cond in conds:
            col = f"detect_{name}_{cs}_{cond}"
            full = res[name][(cs, cond)]["detect"]
            bad = [i for i, r in enumerate(rows) if f"{full[i]:.4f}" != r[col]]
            worst = max(worst, len(bad))
            if bad:
                raise ValueError(f"{col}: {len(bad)} targets do not round to the stored value (first: row {bad[0]})")
    # the summary's medians and shares
    prim = json.loads(SUMMARY.read_text())["primary"]["oracle"]
    for name in ("opm_matched", "opm_dense"):
        for cond in conds:
            r = np.log2(res[name][("opm", cond)]["detect"] / res["squid"][("combined", cond)]["detect"])
            s = prim[f"{name}/combined/{cond}"]
            if abs(float(np.median(r)) - s["median_log2"]) > 1e-9 or abs(float(np.mean(r > 0)) - s["share_opm_better"]) > 1e-12:
                raise ValueError(f"{name}/{cond}: median or share differs from {SUMMARY.name}")

    with open(OUT, "w", newline="") as fh:
        io.csv_status(fh, "DERIVED (full-precision detectability of g2_targets.csv's targets, primary band and noise model "
                          "rebuilt as scripts/g2_band_sensitivity.py does; rounding reproduces the stored table)")
        fh.write("# detect_<array>_<condition>: known-topography detectability of a 10-nAm cortical-normal dipole with the oracle "
                 "covariance (Neuromag: all 306 channels); ratio_<array>_<condition>: OPM over Neuromag; ahead_<array>_"
                 "<condition>: 1 if the ratio exceeds 1 at full precision; values in shortest round-trip decimal form\n")
        wr = csv.writer(fh)
        head = ["hemi", "vertno"]
        for cond in conds:
            head += [f"detect_{ARRAYS[n]}_{cond}" for n in ARRAYS]
            head += [f"{k}_{ARRAYS[n]}_{cond}" for n in ("opm_matched", "opm_dense") for k in ("ratio", "ahead")]
        wr.writerow(head)
        for i, v in enumerate(st.src.target):
            row = [int(st.cortex.hemi[v]), int(st.cortex.vertno[v])]
            for cond in conds:
                ref = res["squid"][("combined", cond)]["detect"][i]
                row += [repr(float(res[n][("combined" if n == "squid" else "opm", cond)]["detect"][i])) for n in ARRAYS]
                for n in ("opm_matched", "opm_dense"):
                    ratio = float(res[n][("opm", cond)]["detect"][i] / ref)
                    row += [repr(ratio), int(ratio > 1)]
            wr.writerow(row)
    n_behind = {cond: int(sum(res["opm_dense"][("opm", cond)]["detect"] <= res["squid"][("combined", cond)]["detect"]))
                for cond in conds}
    print(f"{st.nt} targets; every stored value reproduced; dense array not ahead at {n_behind} targets -> "
          f"{OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
