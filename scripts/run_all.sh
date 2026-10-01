#!/usr/bin/env bash
# Reproduce every milestone in order (see README.md). Requires the MNE sample data in
# data/external/MNE-sample-data. Full-resolution lead fields take ~6 min each; the whole run
# takes roughly 2 h on a laptop (about 40 min of it for the lead fields). Each driver writes its
# outputs to results/<milestone>/ with the code commit and package versions in its JSON.
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
$PY scripts/study_bem_sphere_accuracy.py                 # G2   BEM accuracy near a regular surface (exact sphere test)
$PY scripts/study_opm_near_mesh.py                       # G2   BEM convergence at the OPM cells (head surface 5,120 vs 20,480)
$PY scripts/study_bem_skin_refinement.py                 # G2   headline vs head-surface refinement and number of layers
$PY scripts/g4_epilepsy_adult.py                         # G4   IED detection, adult
$PY scripts/g4_localization.py                           # G4   bounded localization, adult
$PY scripts/study_g4_vs_g2.py                            # G4   detection vs the G2 detectability at the G4 locations
# after the adult baseline (adult-baseline-v1):
$PY scripts/g3a_jas_size_benchmark.py                    # G3A  Jas Table 1 / Fig. 5 (REPRO)
$PY scripts/build_site.py                                # G5   local report in site/_build (not deployed)
