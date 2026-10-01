#!/usr/bin/env bash
# Reproduce every milestone in order (see README.md). Requires the MNE sample data in
# data/external/MNE-sample-data. Full-resolution lead fields take ~6 min each; the whole run
# takes roughly 2 h on a laptop for the adult part (about 40 min of it for the lead fields) and
# about 3 h for the pediatric part. Each driver writes its
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
# after the adult baseline (adult-baseline-v2); needs the 12-, 18- and 24-month infant templates (README, Data setup):
$PY scripts/g3a_jas_size_benchmark.py                    # G3A  Jas Table 1 / Fig. 5 (REPRO)
$PY scripts/g3b_pediatric_helmet.py                      # G3B  fixed adult helmet vs head-adaptive OPM (NEW)
$PY scripts/g4_epilepsy_pediatric.py infant2yr --detection --localization   # G4  IED detection and localization, 2-year template
$PY scripts/g4_epilepsy_pediatric.py infant18mo --detection --localization  # G4  the same, 18-month template
$PY scripts/g4_epilepsy_pediatric.py infant12mo --detection --localization  # G4  the same, 12-month template
$PY scripts/g4_epilepsy_pediatric.py school --detection --localization      # G4  the same, school-age size control
$PY scripts/g4_epilepsy_pediatric.py size2yr --detection --localization     # G4  the same, 2-year size control
$PY scripts/g4_epilepsy_pediatric.py --compare                              # G4  pediatric vs adult comparison and report
$PY scripts/g4_motion.py                                                    # G4  head motion and OPM slippage (bounded extension)
$PY scripts/build_site.py                                # G5   local report in site/_build (not deployed)
