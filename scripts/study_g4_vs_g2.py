#!/usr/bin/env python3
"""Consistency of the G4 detection study with G2 (docs/methods.md section 9).

At the 72 G4 event locations: the median G2 known-topography detectability ratio dense OPM /
Neuromag combined (intrinsic+brain+env, the noise of the time-domain simulation) per G4 depth band,
and the oracle strength-for-50 %-detection ratio Neuromag / dense OPM of the same band. Reads
results/g2/g2_targets.csv and results/g4/g4_adult_summary.json; writes
results/g4/g4_g2_consistency.json.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from opmsquid import io  # noqa: E402

RES = ROOT / "results"
COLS = ("detect_opm_dense_opm_intrinsic+brain+env", "detect_squid_combined_intrinsic+brain+env")


def main():
    g4 = json.loads((RES / "g4" / "g4_adult_summary.json").read_text())
    with open(RES / "g2" / "g2_targets.csv") as fh:
        rows = {(int(r["hemi"]), int(r["vertno"])): r for r in csv.DictReader(fh)}
    by_band = {}
    for loc in g4["locations"]:
        r = rows[(loc["hemi"], loc["vertex"])]
        by_band.setdefault(loc["stratum"][0], []).append(float(r[COLS[0]]) / float(r[COLS[1]]))
    s50 = {k: g4["detectors"][k]["strength_for_50pct_nAm"] for k in ("squid/combined", "opm_dense/opm")}
    out = dict(status="NEW check: G4 detection vs G2 detectability at the G4 locations", inputs=dict(
        g2_targets="results/g2/g2_targets.csv", g4_summary="results/g4/g4_adult_summary.json",
        g4_simulated_at_commit=g4.get("simulated_at_commit")), bands={})
    for b in sorted(by_band):
        sq, op = s50["squid/combined"][f"oracle/depth{b}"]["value"], s50["opm_dense/opm"][f"oracle/depth{b}"]["value"]
        out["bands"][f"depth{b}"] = dict(n_locations=len(by_band[b]), g2_median_ratio_dense_over_combined=float(np.median(by_band[b])),
                                         oracle_s50_ratio_combined_over_dense=None if None in (sq, op) else float(sq / op))
        print(f"depth band {b}: {len(by_band[b])} locations, G2 median ratio {np.median(by_band[b]):.2f}, oracle S50 ratio "
              f"{'-' if None in (sq, op) else f'{sq / op:.2f}'}")
    io.write_json(out, RES / "g4" / "g4_g2_consistency.json")


if __name__ == "__main__":
    main()
