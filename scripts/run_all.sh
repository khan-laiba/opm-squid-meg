#!/usr/bin/env bash
# Reproduce every milestone in order (see README.md). Requires the MNE sample data in
# data/external/MNE-sample-data. Full-resolution lead fields take ~3-10 min each; the whole run
# takes roughly 1.5-2 h on a laptop. Each driver writes its outputs to results/<milestone>/ with the
# code commit and package versions in its JSON.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python

$PY -m unittest discover -s tests -t .                  # physics, units, metrics, geometry, noise
$PY scripts/g1a_jas_benchmark.py                         # G1A  Jas et al. 2026 analytical benchmark (REPRO)
$PY scripts/compute_fullres_forwards.py                  # full-resolution lead fields (all jobs; skips existing)
$PY scripts/g1b_hunold.py                                # G1B  Hunold et al. 2016 (ADAPT)
$PY scripts/g1c_goldenholz.py                            # G1C  Goldenholz et al. 2009 (ADAPT)
$PY scripts/g2_adult_comparison.py                       # G2   realistic adult OPM vs Neuromag (NEW)
$PY scripts/g2_band_sensitivity.py                       # G2   frequency band and OPM response
$PY scripts/g2_report.py                                 # G2   per-configuration report
$PY scripts/g4_epilepsy_adult.py                         # G4   IED detection, adult
$PY scripts/g4_localization.py                           # G4   bounded localization, adult
# after the adult baseline (adult-baseline-v1):
$PY scripts/g3a_jas_size_benchmark.py                    # G3A  Jas Table 1 / Fig. 5 (REPRO)
