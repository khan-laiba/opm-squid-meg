#!/usr/bin/env bash
# Reproduce every milestone in order (see README.md). Requires the MNE sample data in
# data/external/MNE-sample-data. Full-resolution lead fields take ~6 min each; the whole run
# takes roughly 2 h on a laptop for the adult part (about 40 min of it for the lead fields) and
# about 9 h for the pediatric part (G3B, eight pediatric G4 runs of about 1 h each, the motion
# extension). Each driver writes its outputs to results/<milestone>/ with the code commit and
# package versions in its JSON.
#
# Resume: a step that completes leaves a stamp (cache/run_all/<step>.done) holding the code
# commit it ran at. `RESUME=1 scripts/run_all.sh` skips every step stamped at the current commit
# and runs the rest in order, so an interrupted run continues where it stopped; a step whose
# commit differs (code changed since) runs again. Nothing is skipped or stamped on a dirty tree;
# the commit and the tree state are read again before every step.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python
STAMPS=cache/run_all
mkdir -p "$STAMPS"

step() {  # step <name> <command...>
  local name=$1
  shift
  COMMIT=$(git rev-parse --short HEAD)
  DIRTY=$(git status --porcelain -- src scripts configs tests requirements.txt 'legacy/*.py' | wc -l | tr -d ' ')
  if [ "${RESUME:-0}" = 1 ] && [ "$DIRTY" = 0 ] && [ -f "$STAMPS/$name.done" ] && [ "$(cat "$STAMPS/$name.done")" = "$COMMIT" ]; then
    echo "[run_all] $name: completed at $COMMIT, skipped"
    return 0
  fi
  echo "[run_all] $name"
  "$@"
  if [ "$DIRTY" = 0 ]; then echo "$COMMIT" > "$STAMPS/$name.done"; fi
}

step tests $PY -m unittest discover -s tests -t .                       # physics, units, metrics, geometry, noise
step g1a $PY scripts/g1a_jas_benchmark.py                                # G1A  Jas et al. 2026 analytical benchmark (REPRO)
step fullres $PY scripts/compute_fullres_forwards.py                     # full-resolution lead fields (all jobs; skips existing)
step g1b $PY scripts/g1b_hunold.py                                       # G1B  Hunold et al. 2016 (ADAPT)
step g1c $PY scripts/g1c_goldenholz.py                                   # G1C  Goldenholz et al. 2009 (ADAPT)
step g2 $PY scripts/g2_adult_comparison.py                               # G2   realistic adult OPM vs Neuromag (NEW)
step g2_bands $PY scripts/g2_band_sensitivity.py                         # G2   frequency band and OPM response
step bem_sphere $PY scripts/study_bem_sphere_accuracy.py                 # G2   BEM accuracy near a regular surface (exact sphere test)
step near_mesh $PY scripts/study_opm_near_mesh.py                        # G2   BEM convergence at the OPM cells (head surface 5,120 vs 20,480)
step bem_skin $PY scripts/study_bem_skin_refinement.py                   # G2   headline vs head-surface refinement and number of layers
step head_surface $PY scripts/study_head_surface_effect.py               # G2   headline under the v3 and v4 head models (A-BEM-CONFORM)
step g2_report $PY scripts/g2_report.py                                  # G2   per-configuration report (reads the studies above)
step g4_adult $PY scripts/g4_epilepsy_adult.py                           # G4   IED detection, adult
step g4_loc $PY scripts/g4_localization.py                               # G4   bounded localization, adult
step g4_vs_g2 $PY scripts/study_g4_vs_g2.py                              # G4   detection vs the G2 detectability at the G4 locations
# after the adult baseline; needs the 12-, 18- and 24-month infant templates and the school-aged
# children (scripts/fetch_school_subjects.py; README, Reproduce):
step g3a $PY scripts/g3a_jas_size_benchmark.py                           # G3A  Jas Table 1 / Fig. 5 (REPRO)
step school_prep $PY scripts/prepare_school_subjects.py                  # G3B  school-aged children: modelled skull, fiducials, source space
step child_bem $PY scripts/study_child_bem.py                            # G3B  the modelled skull (A-BEM-CHILD) checked on the adult
step school_checks $PY scripts/study_school_anatomy.py                  # G3B  skull depth, fiducial transfer and talairach checks
step g3b $PY scripts/g3b_pediatric_helmet.py                             # G3B  fixed adult helmet vs head-adaptive OPM (NEW)
step g4_infant2yr $PY scripts/g4_epilepsy_pediatric.py infant2yr --detection --localization    # G4  IED detection and localization, 2-year template
step g4_infant18mo $PY scripts/g4_epilepsy_pediatric.py infant18mo --detection --localization  # G4  the same, 18-month template
step g4_infant12mo $PY scripts/g4_epilepsy_pediatric.py infant12mo --detection --localization  # G4  the same, 12-month template
step g4_school $PY scripts/g4_epilepsy_pediatric.py school --detection --localization          # G4  the same, school-age size control
step g4_size2yr $PY scripts/g4_epilepsy_pediatric.py size2yr --detection --localization        # G4  the same, 2-year size control
step g4_childA $PY scripts/g4_epilepsy_pediatric.py childA --detection --localization          # G4  the same, school-aged child A
step g4_childB $PY scripts/g4_epilepsy_pediatric.py childB --detection --localization          # G4  the same, school-aged child B
step g4_childC $PY scripts/g4_epilepsy_pediatric.py childC --detection --localization          # G4  the same, school-aged child C
step g4_compare $PY scripts/g4_epilepsy_pediatric.py --compare           # G4  pediatric vs adult comparison and report
step g4_matched $PY scripts/study_g4_matched_rate.py                     # G4  detection at matched held-out false-event rates
step fit_failures $PY scripts/study_g4_fit_failures.py                   # G4  failed fits under declared criteria (from the event tables)
step g4_motion $PY scripts/g4_motion.py                                  # G4  head motion and OPM slippage (bounded extension)
step site $PY scripts/build_site.py                                      # G5   local report in site/_build (not deployed)
