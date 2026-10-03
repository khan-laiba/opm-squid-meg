# PLAN — adult-to-pediatric OPM vs SQUID MEG study

Milestone plan for `GOAL.md`. Each task is labelled:

- **REPRO**: reproduction of a published result, with the paper's own definitions.
- **ADAPT**: methodological adaptation (original anatomy, recordings or details unavailable).
- **NEW**: new experiment or study choice.

Evidence (commit, command, outputs) is logged per milestone. Parameters and their provenance
live in `docs/provenance_register.md`; methods in `docs/methods.md`.

## Status

| Milestone | Status | Evidence |
|---|---|---|
| G0 audit, provenance, plan | done; private repository `khan-laiba/opm-squid-meg`, tag `g0` | this file; `docs/audit.md`; `docs/literature/`; `docs/provenance_register.md` |
| G1A Jas analytical benchmark | done; internally reviewed (approve with notes; notes addressed) | `scripts/g1a_jas_benchmark.py` -> `results/g1a/`; Eq. 1 vs Sarvas 2-D maximum 4.7e-15, vs MNE sphere 5.4e-8; d_eq(eta 3) = 27.665 mm |
| G1B Hunold depth-orientation spikes | done; internally reviewed (approve with notes; notes addressed); rerun with the 4-mm source rule and a realization-averaged Fig. 6 calibration | `scripts/g1b_hunold.py` -> `results/g1b/`; calibrated background (scalar 0.47; one realization alone 0.45-0.52), p2p: bin means 0.74-0.85x the paper (0.64-0.94x over the calibration range), noisy p2p 0.94-0.99x, r 0.93-0.98, 2.5-classification agreement 79-89 %, GM-MM sign agreement 53/54 dipole and 49/49 patch bins; v4 (computed at e53bea8): sparse patch bins (< 5 patches) left blank, segment-only Hilbert variant reported |
| G1C Goldenholz cortical SNR maps | done; internally reviewed (approve with notes; notes addressed); rerun with the 4-mm source rule | `scripts/g1c_goldenholz.py` -> `results/g1c/` (computed at e53bea8); s_s 1.77 nAm at 4,000 sources (paper 1.6-1.9); focal median -22.1 dB, 56 % inside the paper's -29/-19 dB range (skull 0.006 S/m, the probable intended value; 0.06 S/m as printed run and labelled); deep medial cortex darkest; without the medial wall and with untruncated patches (v4 variant) the medians rise by 0.5-0.8 dB |
| G2 realistic adult OPM-Neuromag comparison | done; internally reviewed; `adult-baseline-v1` (1.21x) corrected in `adult-baseline-v2` (near-surface BEM error), `adult-baseline-v3` (whole-cell OPM clearance) and `adult-baseline-v4` (every anatomy's head surface on its MRI scalp: equal OPM standoff; see the decisions log) | `scripts/g2_adult_comparison.py` -> `results/g2/` (`G2_report.md`; computed at e53bea8): with modelled brain noise the dense 208-site OPM array is 1.14x [1.12-1.17] Neuromag combined (1.01-1.30x for OPM noise 30-7 fT/sqrt(Hz), 0.92-1.14x jointly with a 0-6 mm scalp gap; 3- and 1-layer BEM within 1.3 %) and 1.12x [1.08-1.15] after the 8-term external-field projection (v3: 1.05x, a tie; the change comes from deep sources, `head_surface_effect.json`), the matched 98-site array 1.01x (a tie; 0.95x after the projection), a triaxial OPM at the matched sites (294 channels) 1.12x (1.05x with doubled tangential noise); 1.4-1.6x for sources within 20 mm of the scalp, 1.05-1.07x below 35 mm; peak-channel SNR favours Neuromag (0.90x); without brain noise Neuromag wins (dense 0.74x); every head-model- and array-dependent adult result rerun at ed852b2 (detection, BEM studies) or e53bea8 (G1B, G1C, G2, localization), the exact-sphere BEM check (no array, no head model) at 81168f3 |
| G3 pediatric extension | G3A done (REPRO size benchmark); G3B done after the adult freeze (goal review 2026-10-02: partially met, see `docs/goal_review.md`; its adult-only OPM standoff removed in v4, where every head surface lies on its MRI scalp (A-BEM-CONFORM); native school-aged anatomy added 2026-10-03, with a modelled skull) (24-, 18- and 12-month infant templates, three school-aged children of OpenNeuro ds005234, adult scaled to school-age and 2-year head size; fixed Neuromag helmet at source-blind placements, yaw included, vs refitted OPM; counterfactual helmet; 5- to 20-mm and fixed-density patches) | `scripts/g3a_jas_size_benchmark.py` -> `results/g3a/`; `scripts/g3b_pediatric_helmet.py` -> `results/g3b/` (`G3B_report.md`; computed at e53bea8, every array at a median 7.0 mm above its MRI scalp): dense OPM vs Neuromag combined at top contact, D_adult +1.00 dB, Delta +0.44 (school-age size), +1.07 (2-year size), +0.73 (2-year template), +0.88 (18 months), +0.96 dB (12 months) (v3, with the adult's sensors 0.8 mm farther out: +0.59 to +1.25 dB); with a helmet scaled with the head about the laterally centred head -0.21, -0.38, -0.30, -0.16, -0.18 dB (the gain is the fixed helmet's fit); robust to OPM noise 7-30 fT/sqrt(Hz), background x0.5/x2, 1-layer BEM (the school-aged children's difference of the medians turns slightly negative at 7 fT/sqrt(Hz), -0.15 to -0.03 dB); top contact is the adult's second-lowest of 12 source-blind placements (against each head's family median the differences of the medians are 0.14-0.20 dB smaller); school-aged children A-C (7.8-8.7 years, head circumference 520, 486, 535 mm; computed at 71de176, the six earlier anatomies reproduced exactly): Delta +0.30, +0.41, +0.17 dB (depth-reweighted -0.04, +0.05, -0.05 dB; projected +0.03, +0.26, -0.19 dB), counterfactual about the laterally centred head -0.76, -0.59, +0.04 dB: at a similar head circumference children A and C gain less than the scaled adult (+0.44 dB; intervals overlapping or touching), child B (smaller head) about as much; the redraws at 6ebac95 and ff6e25c mark child C's infeasible x-5mm placement and correct the report's notes, every number unchanged |
| G4 epilepsy detection and localization | adult and pediatric (three templates, both size controls, three school-aged children) done; v4: detection rerun at ed852b2, localization at e53bea8 (with MNE's own dSPM, Neuromag magnetometers or gradiometers alone and a held-out check of the thresholds), the cross-reading summaries at 9512f1b; the school-aged children (both studies) at 71de176, the summaries again at 011cec1 (earlier anatomies unchanged); motion and slippage bounded extension done (ed852b2) | `scripts/g4_epilepsy_adult.py`, `scripts/g4_localization.py` -> `results/g4/`; 50 % detection at 1 false event/min: dense OPM 34 vs Neuromag combined 47 nAm at 10-20 mm (paired ratio 1.36 [1.10-1.58]; 11/1 locations, p = 0.004; oracle 15/0, p = 0.00006; v3 8/2, ratio 1.29; v2 16/0, 1.51), 254 vs 302 nAm at 45-70 mm; v4 also favours the dense array at 20-30 mm (9/1) and 45-70 mm (9/0), uncorrected and found in only one earlier run (20-30 mm, an early v2 run): reported, not established; matched: not ahead in any band (7/4 at 10-20 mm, p = 0.40) and no longer behind at 45-70 mm (3/6; v3 1/9, v2 0/7); localization: ECD similar (~4-6 mm); the dSPM errors for 320-nAm patches are 5.6 and 2.4 mm smaller (p = 0.054 and 0.067; with MNE's dSPM 1.0 and 2.3 mm; v3 0 mm, v2 ~5 mm): not robust; pediatric (`scripts/g4_epilepsy_pediatric.py`, `G4_pediatric_report.md`): 10-20 mm strength ratio Neuromag/OPM 1.33-1.65 in the templates and size controls, 1.18-1.43 in the school-aged children (adult 1.36; their p = 0.004-0.022 does not survive the correction over 24 comparisons), deeper differences only for the 2-year size control and the 12- and 18-month templates among the children (uncorrected; none survives a correction); localization comparisons with p < 0.05 favour an OPM array 64 to 1 (of 360 in nine anatomies; uncorrected, about 18 expected by chance, 28 for weak 80-nAm sources; none of the school-aged children's survives a correction); motion (`scripts/g4_motion.py`, `G4_motion_report.md`): uncompensated head displacement costs Neuromag 1.4-1.7 dB at 10 mm, a 3-deg OPM cap slip 0.1-0.5 dB; un-modelled in-band rotation costs the OPM 1 dB at ~0.02 deg RMS in a 1-nT field (no correction) or 0.4-0.5 deg (field projection, 1 deg/1 % calibration), thresholds inversely proportional to the field |
| G5 repository, tests, report | private repository `khan-laiba/opm-squid-meg` (no Pages, verified 2026-10-01); 160 unittest tests pass with the v4 results and the school-aged children (a clean clone of 8bdeb8b passed 159 of its 160, one skip without the local lead-field cache; db3f0dc 147 of 148; v3: a clean clone of 9771f68 passed 132 of its 133); `scripts/run_all.sh`; local report (`scripts/build_site.py` -> `site/_build/`, link-checked, not deployed); clean-environment smoke test passed (2026-09-30, at a commit before 81168f3, 87 tests then); release prepared, not executed | `docs/release_checklist.md`: figures in open fonts, no local paths in output summaries; owner decisions left: licence, git author metadata, template/fsaverage-derived figures, the CC-BY raster, two documents naming local folders, visibility and Pages |

## Decisions log

| Date | Decision | Source |
|---|---|---|
| 2026-09-30 | `GOAL.md` = verbatim copy of `OPM_SQUID_Adult_to_Pediatric_Goal.md`; condensed `/goal` text in `goal_condition.txt` (3,352 characters). | owner |
| 2026-09-30 | New PRIVATE GitHub repository `khan-laiba/opm-squid-meg` (none existed); push after G0; no Pages deployment or public visibility without explicit approval. | owner |
| 2026-09-30 | Real Neuromag geometry from the MNE sample dataset (`MNE-sample-data-processed.tar.gz`, 1.58 GB, osf.io via `mne.datasets.sample`), stored in `data/external/` (not committed). | owner |
| 2026-09-30 | Earlier scripts/outputs moved unchanged to `legacy/`; the conformal dense arrays are labelled an idealized baseline. | GOAL G0 |
| 2026-09-30 | Primary adult anatomy for G1B/G1C/G2: MNE `sample` subject (individual adult MRI, BEM surfaces, real head position in a Vectorview helmet). fsaverage is a secondary template. | NEW (see Inputs) |
| 2026-09-30 | G1B: p2p numerator primary (Fig. 6 digitisation); two background levels (as specified; one scalar from the Fig. 6 magnetometer baselines, gradiometers as independent check); paper maps never used for calibration. | NEW (U-HU-numerator, U-HU-bglevel) |
| 2026-09-30 | Channels marked bad in the sample recording (MEG 2443) are excluded wherever recorded noise enters (G1C; G2 measured-noise calibration). | NEW |
| 2026-09-30 | Pediatric anatomy for G3B: the 2-year infant template (`mne.datasets.fetch_infant_template('2yr')`, `ANTS2-0Years3T.zip`, 391,756,700 bytes, github.com/christian-oreilly/infant_template_paper, LGPL-2.1), plus the adult scaled to school-age and 2-year head size as size-only controls. No school-aged native anatomy (another download would be needed). | owner (download approved) |
| 2026-09-30 | Clean-environment smoke test: fresh venv from `requirements.txt` (pinned packages downloaded with approval), clone of the repository, unit tests and G1A: passed (87 tests; G1A identical apart from provenance). | owner (download approved) |
| 2026-09-30 | `adult-baseline-v2` results recomputed at one clean commit (81168f3) after the second independent review (cell clearance checked in the forward model's frame). | NEW (review v2) |
| 2026-09-30 | `adult-baseline-v2` tagged at c472694 after a light independent re-review (approve with notes, fixed) and pushed; repository private, no Pages (verified). | owner rules |
| 2026-09-30 | G3B design: primary placement = top contact (20 mm), centred/back/bounded variants as sensitivity; counterfactual helmet scaled with the head-circumference ratio about the head origin; D/Delta in dB of known-topography detectability; vertex-wise homology for scaled controls, parcels and strata for the template; medial wall excluded from G3B summaries. Developed on branch `g3b` (worktree) while the adult reruns ran; merged with --no-ff. | NEW (A-G3-*) |
| 2026-10-01 | More pediatric anatomies: the 18- and 12-month templates of the same series (`fetch_infant_template('18mo' / '12mo')`, `ANTS18-0Months3T.zip` 385,415,605 bytes and `ANTS12-0Months3T.zip` 377,092,022 bytes, same repository and licence) in G3B and pediatric G4; pediatric G4 also for the 2-year size control. | owner (downloads approved) |
| 2026-10-01 | G4 motion extension (bounded, secondary): sustained displacement in the fixed helmet vs OPM cap slip (geometry known vs mismatched template), in-band room-field coupling of the head-mounted array (no correction, homogeneous or 8-term projection; calibration errors; artefact outside and inside the noise model), exact rigid motion over a recording; per unit field. | NEW (A-MOT-*) |
| 2026-10-01 | Release preparation, staying private: figures in open fonts by default (licensed artwork fonts only on request, for the pixel verification), output summaries without local absolute paths, `docs/release_checklist.md`; no visibility or Pages change. Provenance now also watches `legacy/*.py` (G1A imports the Fig. 3 replica). | owner ("prepare, stay private") |
| 2026-10-01 | Final end-to-end reviews of every finding by two other models (owner request): Fable 5.1 (approve with notes; addressed in 90c18ae and 55a116d) and GPT 6 Astra via Codex (reject pending revision: whole-cell OPM clearance, detection scoring, patch support, censored intervals, report-builder deletion, CSD phase, ROC count, fingerprints). All findings fixed (section "Final reviews" below) and every OPM-dependent result recomputed at one clean commit (`adult-baseline-v3`). | owner |
| 2026-10-02 | End-to-end goal review of G0-G5 (`docs/goal_review.md`; three independent reviewer agents, a clean-clone test run): G0, G1, G4 and G5 met for the private deliverable, G2 met with reservations, G3 partially met (the clearance rule left the adult's OPM sensors about 0.8 mm farther from the scalp than the children's: equalised in a sensitivity (`scripts/study_g3b_standoff.py`), the fixed-helmet Delta stays +0.47 to +1.12 dB and the counterfactual Delta is -0.34 to -0.06 dB, the same reading; no school-aged anatomy, owner-blocked). Fixed in the review: failed fits counted, the two G1C side-effects measured, configurations as the source of truth, the metric dependence and the deep reversal after the projection in the reports, documentation versions, the smoke command, disclosures. | owner |
| 2026-10-01 | Re-check of the revision (0d6b992): Codex (reject pending focused corrections: censored S50-ratio intervals, the clearance guarantee, MNE's standard coil definitions in the fingerprints) and an independent verifier agent (approve with fixes: ten text corrections). All fixed; the G4 summaries re-summarised from the stored simulations (d667a9d), the lead fields recomputed. | owner |
| 2026-10-02 | Equal OPM standoff (owner request after the goal review): every anatomy's BEM head surface has its vertices on its MRI scalp (A-BEM-CONFORM; the templates' are built so, the sample's stored outer skin lay a median 0.8 mm outside its scalp over the sensor region). This removes two adult-only differences in G3B's D_adult: the clearance rule's outward shift (v3 median sensor height 7.78 / 7.76 mm vs 7.00 mm in every child; v4 6.99-7.01 mm everywhere) and the adult's axis tilt from its outer skin's stored normals (A-OPM-AXIS). The rule is defined by geometry alone (the templates' construction applied to the adult) and was fixed before any v4 result existed; it changes the adult baseline after the pediatric outcomes were known, a deviation from the freeze sequence stated here and in the goal review addendum. Every result that depends on the head model or the OPM arrays recomputed (`adult-baseline-v4`); the standoff sensitivity study is retired. | owner |
| 2026-10-02 | The pending review items (goal review "open"): G1B minimum-count rule for patch bins and the segment-only Hilbert variant; G1C variant without the medial wall and with untruncated patches; G2 plug-in covariances from one common noise realization, Neuromag noise from the measured empty-room spectrum (A-G2-SQUIDMEAS) and a triaxial channel-count control (A-OPM-TRIAX); G3B yaw placement variant, 20-mm and fixed-density patches (A-G3-PATCH) and the placement citations (A-G3-PLACE); G4 localization with Neuromag magnetometers or gradiometers alone, dSPM through `mne.minimum_norm` next to the study's own, and the detector thresholds on held-out null data; status lines in every result CSV and on the G1A replica's summary; `RESUME=1 scripts/run_all.sh`. Developed on branch `pending-items` (worktree) while the equal-standoff rerun ran; the G1B/G1C/G2/G3B/localization results come from its commit. Not done (owner decisions): a CI workflow; a school-aged native anatomy (needs an approved download). | owner |
| 2026-10-03 | School-aged native anatomy (owner approvals): three typically developing children of OpenNeuro ds005234 v2.2.0 (Fadeev et al. 2024; sub-Z213, sub-Z209, sub-Z226; 7.8-8.7 years): 39 files, 122,745,464 bytes (white, sphere and dense scalp surfaces, aparc annotations, and the talairach.xfm and BEM surfaces stored in their folders), then 6 files, 2,215,848 bytes (children B's and C's own watershed BEMs; the snapshot's file tree shifts each subject's BEM and talairach.xfm into the preceding subject's folder, so child A's came with the first files, in sub-Z209's folder). The watershed inner skull lies just below the scalp, so the skull is modelled (A-BEM-CHILD); no fiducials come with the data, and fiducials from MNI coordinates miss the adult's own by 10-27 mm (and only child A's own talairach.xfm was obtained), so the adult's are transferred (A-G3-FID); both checked (on the adult and on the templates: `scripts/study_child_bem.py`, `scripts/study_school_anatomy.py`). Added to G3B and the pediatric G4 studies as children A-C; the earlier anatomies' results reproduce unchanged. | owner |

## Inputs and availability

| Input | Needed for | Status |
|---|---|---|
| Jas et al. 2026 preprint (PDF) | G1A, G3 size benchmark | available locally (not committed) |
| Hunold et al. 2016 (PDF) | G1B | available locally (not committed) |
| Goldenholz et al. 2009 (PDF) | G1C | available locally (not committed) |
| TRIUX specification image | G2 noise and geometry references | **not found on disk** (project folders searched again 2026-10-01; personal folders not searched); values transcribed in GOAL.md are used and flagged |
| MRN Neuromag page (T3 coils 3014/3024, shielded-room floor 5-7 fT) | G2 | read 2026-09-30 |
| MNE implementation docs (coil definitions, frames), MNE 1.13.2 | all | read 2026-09-30 |
| Neuromag 3-D geometry, `dev_head_t`, head-MRI trans | G2 | MNE sample data (MGH Vectorview: 204 grads coil 3012 T1, 102 mags coil 3024 T3) |
| Measured SQUID noise | G2 measured-noise scenario | MNE sample `ernoise_raw.fif` (empty room) and baseline covariance |
| Adult anatomy | G1B, G1C, G2 | MNE `sample` subject (3-layer BEM surfaces, oct-6 source space); fsaverage (sibling folder, read-only) |
| School-aged anatomy (FreeSurfer surfaces, BEM surfaces, scalp) | G3 | **obtained 2026-10-03** (route (b): three children of OpenNeuro ds005234, `data/external/school_subjects/`, not committed; skull modelled, fiducials transferred, see the decisions log). Earlier (owner decision 2026-09-30: route (c) + (d)): FreeSurfer is not installed; MNE only packages infant templates (up to 2 years, Neurodevelopmental MRI Database via `mne.datasets.fetch_infant_template`). Routes: (a) age-specific average templates of the Neurodevelopmental MRI Database (registration required); (b) an OpenNeuro school-aged dataset that ships FreeSurfer outputs (e.g. fMRIPrep `sourcedata/freesurfer`); (c) the 2-year infant template as a young-child anatomy (**downloaded**, `data/external/infant_subjects/ANTS2-0Years3T`), and since 2026-10-01 the 18- and 12-month templates of the same series (**downloaded**); (d) scaled adult surfaces as a labelled size-only control (**used**: school-age and 2-year size). Between-child variability is shown by three individuals of one dataset, not estimated. |
| Original Hunold/Goldenholz participant data and recordings | exact reproduction | unavailable; G1B and G1C are adaptations |
| Jas SEF recordings | experimental validation | unavailable; used as context only |

## G0 — audit, provenance, plan

- [x] Goal files, condensed goal, git repository, pinned `requirements.txt`.
- [x] Legacy work moved to `legacy/` with a status table (`legacy/README.md`).
- [x] Full-text extractions: `docs/literature/{jas2026,hunold2016,goldenholz2009}.md` (each independently verified).
- [x] `docs/audit.md`: validated vs unverified legacy claims; reusable components.
- [x] `docs/provenance_register.md`: paper-reported (J, HU, GO), hardware (HW), recovered (R), new assumptions (A), decisions (D); unresolved ambiguities listed.
- [x] `docs/methods.md` (living document).
- [x] Private repository created, first push, tag `g0`.

## G1 — adult foundations

### G1A Jas analytical benchmark (REPRO)

1. `opmsquid.sphere`: Eq. 1 peak radial field; independent Sarvas/Biot-Savart numerical field
   with a 2-D maximisation over the whole sensor sphere; MNE sphere forward as a third check.
2. Fig. 3 (reuse the validated legacy replica) and Fig. 4 (eta = 1-6: SNR curves and d_eq vs eta).
3. Exclusions: r_Q = 0 (sphere centre, silent) and any zero-field source removed from ratios.
4. Absolute noise: document sigma_SQUID = 0.3546 pT (recovered, not stated) and show that
   eta alone fixes only ratios.
5. Toy target-vs-background depth experiment as a separate explanatory benchmark; list every
   text/caption inconsistency with both readings.
6. Checks: analytic vs numerical maxima (relative error target < 1e-6); d_eq(eta = 3) = 27.665 mm
   (paper: about 28 mm); d_eq(eta) existence range (113/95)^3 < eta < 5.3086.

### G1B Hunold depth-orientation spike simulations (ADAPT; OPM part NEW)

Done (`scripts/g1b_hunold.py`, methods section 6, `configs/hunold_reference.toml`):
1. [x] Sample subject, full-resolution white surface, 3-layer BEM with the paper's conductivities,
   4-point coil integration; depth/orientation from BEM nodes; the paper's bins; 3783 dipoles
   stratified to the paper's per-bin counts; 20-mm^2 patches (2928 grow).
2. [x] Background: 10 % of vertices, band-limited, stationary (edge-transient defect fixed),
   +/-10 nAm per dipole; one realization shared by sources and arrays.
3. [x] SNR as defined by the paper with the numerator inferred from Fig. 6 (the paper prints "spike
   peak" without saying noise-free or noisy, peak or peak-to-peak; U-HU-numerator): p2p primary,
   peak, noisy peak and noisy p2p as variants.
4. [x] Absolute level: the paper's Fig. 6 baselines are 0.47x (MM) / 0.42x (GM) our expected
   baselines (20 realizations); both levels reported (U-HU-bglevel). Calibrated, p2p: strong
   bins 0.87-0.92x the paper, weak bins 0.73-0.90x; noisy p2p 0.97-1.04x (the paper lies between
   the two numerators). GM > MM superficially, convergence with depth: reproduced.
5. [x] OPM (matched 95 sites, NEW; v3 arrays): brain noise only, OPM > MM for superficial sources
   (1.4x at 20-25 mm), equal near 40-45 mm, <= MM deeper (0.9x), < GM everywhere (0.69-0.80x, best
   single channel); intrinsic sensor noise (3.5-30 fT/sqrt(Hz)) is negligible against this
   background.
6. [x] Independent review (approve with notes). Fixes: stationary background, calibration range and
   effective-level wording, noisy-p2p variant, spike comparison, usable sources (A-BEM-DIST).

### G1C Goldenholz cortical SNR maps (ADAPT)

Done (`scripts/g1c_goldenholz.py`, methods section 7):
1. [x] 10-nAm dipoles at all valid vertices; 10/16-mm geodesic patches at 50 pAm/mm^2 on the oct-6
   centroids; signed sums.
2. [x] Eq. 1 with 1/N per channel set (mag 102, grad 203, pooled 305; MEG 2443 excluded);
   modelled noise (paper's calibration rule) and recorded noise (task baselines, ADAPT).
3. [x] Skull 0.006 and 0.06 S/m both run and labelled.
4. [x] Soft comparisons with the paper (source SD, display range, medial maps, patch-size effect).
5. [x] OPM extension (NEW), brain-only calibration and OPM noise sweep after review.
6. [x] Independent review (approve with notes). Fixes: near-skull BEM artefacts (A-BEM-DIST),
   double-counted instrument noise in the extension, wording, tests.

GATE G1: all three run reproducibly; each output states REPRO/ADAPT status and assumptions.

## G2 — realistic adult OPM vs Neuromag (NEW)

1. Sensors: sample-data Neuromag geometry with the real `dev_head_t` and head-MRI transform;
   coil types set to MRN T3 (3014 grads, 3024 mags; the 3012-to-3014 change is documented).
   Accurate coil integration. Compute actual scalp-to-pickup-coil gaps per sensor.
2. OPM array: 102 matched sites (Neuromag sensor sites projected to the scalp), sensing
   centre 7 mm from the helmet's inner surface (Jas p. 10, assumption), 10-mm cubic cell
   via a custom coil definition (finite-volume integration tested), extra scalp gap as a
   separate parameter; packing/clearance check. Channel-budget and full-system controls.
3. Noise: intrinsic white noise (mags 3.5 fT/sqrt(Hz), grads 3.6 fT/cm/sqrt(Hz), OPM sweep
   7-30 fT/sqrt(Hz)); cortical background with area-scaled variance (independent baseline
   plus bounded spatially correlated extension), calibrated against the MNE sample
   brain-noise covariance (baseline minus empty-room) as a labelled option; environmental
   interference via an external multipole basis fitted to the empty-room recording.
   PSD-to-variance through the actual composite filter.
4. Metrics: paper-specific (Hunold, Goldenholz) plus peak-channel SNR, mean-power SNR_dB,
   known-topography detectability with rank-aware whitening; oracle vs estimated covariance.
5. Endpoints: signal vs depth, depth-orientation heatmaps, cortical SNR maps, patch size,
   regional relative-performance maps, helmet fit, source-blind head-position sensitivity,
   convergence (sensor integration, BEM, source space) and uncertainty.
6. [x] Freeze the adult baseline as git tag `adult-baseline-v1` (after the independent review's
   pre-tag items P1-P4 and the reruns on the final arrays).

## G3 — pediatric extension (NEW; size benchmark REPRO)

1. [x] REPRO Jas Table 1 / Fig. 5 size-following benchmark (s = h + 18 mm), plus an idealised
   concentric fixed-shell contrast (NEW).
2. [x] Fixed adult Neuromag helmet vs OPM refit per head (`opmsquid.pediatric`,
   `scripts/g3b_pediatric_helmet.py`): the 2-year infant template in native dimensions, and the
   adult scaled to school-age (85/95) and 2-year head circumference as labelled size-only
   controls; since 2026-10-03 three school-aged children with individual MRIs (OpenNeuro ds005234;
   modelled skull, transferred fiducials). Three individuals of one dataset and average templates of
   one database support no claim about anatomical variability.
3. [x] Centred, top- and back-contact and bounded translated/rotated placements chosen
   source-blind; regional gaps (helmet selections), source-to-sensor distances, OPM coverage and
   packing, noise composition; counterfactual helmet scaled with the head (labelled mechanistic
   control).
4. [x] D_child, D_adult, Delta for each SQUID comparator, vertex-wise (scaled controls) and per
   parcel and depth/orientation stratum (templates and school-aged children), area-weighted, parcel-bootstrap intervals,
   sparse strata shown; projection, matched array, patches, secondary metrics; OPM noise,
   background and BEM sensitivity; usefulness maps.
5. [x] Independent review of G3B and pediatric G4 (approve with notes): no code or provenance
   problem; the mechanism statements were wrong (both systems gain in the smaller heads, the OPM
   more; the counterfactual Delta depends on the SQUID comparator, so it is not explained by the
   OPM's fixed standoff alone) and the template's counterfactual residual is confounded by its
   lateral offset. A re-check found six remaining text defects (depth apportionment, the 18-mm
   placement, three numbers, G4 depths), fixed. Fixed: wording, absolute detectability per system, orientation strata, a laterally
   centred placement with its own counterfactual, an 18-mm contact bound, a channel-count control.
6. [x] The refactored G4 scripts (`Context`, 225ffc4) reproduce the adult G4 results of 81168f3
   exactly: the detection and localization summaries are identical apart from provenance and the
   event tables byte-identical (checked 2026-09-30, 21:04 and 21:38; G4-relevant code changed later
   only in `plotting.py`).
7. [x] 18- and 12-month templates added (2026-10-01; G3B computed at cb1b8a9, summaries redrawn at
   fd40dfe): the four anatomies of the reviewed pass reproduce to 1e-9 (the bounds of 188 secondary
   bootstrap intervals moved with the shared random stream); Delta +0.85 [+0.49, +1.27] and +1.03 [+0.66, +1.44]
   dB on the v2 arrays (v3, 19a8fd2: +0.99 [+0.58, +1.33] and +1.23 [+0.73, +1.61] dB; v4, e53bea8: +0.88 [+0.54, +1.05] and +0.96
   [+0.61, +1.51] dB); the depth checks of the earlier hand computation are now computed for every template
   (`template_depth_checks`).

## G4 — epilepsy (NEW)

1. [x] (adult) IED-like events across regions/depths/orientations/extents and strengths.
2. [x] (adult) Oracle vs practical detector; thresholds calibrated on independent null data and
   frozen; sensitivity vs false events/min; paired comparisons with the location as the unit.
3. [x] (adult) Bounded ECD and MNE/dSPM localization with off-grid truth and registration/model
   mismatch (coregistration draws shared by all arrays); nondetections and errors reported
   separately; GOF descriptive only (whitened GOF depends on the channel count).
4. [x] Pediatric runs (`scripts/g4_epilepsy_pediatric.py`): the adult detection and localization studies
   unchanged on the three templates, both size controls and the three school-aged children (Neuromag
   at top contact, refitted OPM arrays); comparison in `results/g4/G4_pediatric_report.md` (regenerated
   from the v4 runs: detection simulated at ed852b2, localization at e53bea8, the school-aged children
   at 71de176, compared at 011cec1; v3 at 19a8fd2/d667a9d, earlier at fd40dfe).
5. [x] Motion and slippage, a bounded secondary extension (`scripts/g4_motion.py`,
   `src/opmsquid/motion.py`, `results/g4/G4_motion_report.md`; methods section 12).
6. [x] End-to-end goal review (2026-10-02, `docs/goal_review.md`): failed fits counted under declared
   criteria (`scripts/study_g4_fit_failures.py`), the detector's and localization's oracle-like
   conveniences and the dSPM implementation's deviation from MNE disclosed, coverage limits stated,
   the pediatric report's held-out range computed from the data, pediatric ROC figures and the
   operating-point notes on the report site.

## G5 — software, reproduction, report

Modular package `src/opmsquid`, scripts per milestone, `configs/`, `tests/` (unittest; the
venv has no pytest), cached forward models keyed by input hashes, fixed seeds, resumable runs,
clean-environment smoke test (needs a package download; ask first), local static report in
`site/` (jinja2), never deployed without approval.

## `adult-baseline-v2`: fixes of the `adult-baseline-v1` review items

1. [x] Sensor-side BEM accuracy. v1's OPM cells reached into the 5,120-triangle head surface
   (integration points up to 2.4 mm inside), where the 3-layer BEM field is not converged on this
   head (an exact sphere test shows the code itself is accurate near a regular surface). v2 refines
   the 3-layer head surface to 20,480 triangles (A-BEM-SKIN) and keeps every integration point
   >= 1 mm outside it (A-OPM-CLEAR). Effect: the dense-array headline fell from 1.21x to 1.13x;
   the 3- and 1-layer models now agree (1.129x vs 1.127x on the convergence subset, final arrays) and
   cell vs point differs by <= 2.8 % in the peak channel (median 0.03 %).
2. [x] Scalp-gap variants move the primary sites outward (same sites) instead of rebuilding.
3. [x] G4 localization: location i uses coregistration draw i mod 8 (all draws in every condition).
4. [x] Provenance: full-resolution lead fields carry fingerprint sidecars (sensors, transforms,
   BEM conductivities and geometry) checked on every load (G1B, G1C, G2, G4); the G4 detection
   summary records its simulation commit.
5. [x] Tests: array composition (97/204/212), cell clearance (exact distances on MNE's own coil
   integration points), MRI-scalp clearance, gap variants, lead-field fingerprints, exact
   point-to-mesh distance.
6. [x] G4 locations are never on the medial wall.
7. [x] Second independent review (approve with notes): the cell-clearance check built each cell in
   the MRI frame while the forward model builds it in the head frame (in-plane axes rotated by up
   to ~100 deg); now checked in the head frame with exact point-to-triangle distances (every
   integration point >= 1.000 mm outside). Arrays: matched 97, dense 212 (was 211), 204-site subset
   of the 212. Every adult result recomputed at one clean commit (81168f3); documentation fixes
   (held-out false-event rates, cell-vs-point wording, coverage counts, the h^2 convergence
   assumption, uncorrected p-values, paired S50 ratios instead of overlapping intervals). The three
   Neuromag full-resolution lead fields (inputs unchanged in v2, so kept by fingerprint) were
   recomputed at 81168f3 without the chunk cache: bit-identical to the matrices the reruns used; every
   lead-field sidecar now records 81168f3.

## Final reviews (2026-10-01): fixes and `adult-baseline-v3`

Codex (GPT 6 Astra) findings, each verified before fixing, and what changed:

1. [x] M1 OPM clearance covered the 27 integration points, not the 10-mm cell: corners reached up to
   0.89 mm inside the BEM head surface (57 dense, 19 matched sites). A-OPM-CLEAR v3 keeps the whole
   cell >= 1 mm outside: integration points and the cell surface sampled at 1-mm spacing (with the
   1-Lipschitz distance the whole surface stays >= 0.29 mm off the head surface; a 41 x 41 check of
   the final arrays gives >= 1.002 mm). Arrays: matched 95 of 102 sites (was 97), dense 205 (was
   212), the 204-site control a subset of the 205; 76 matched and 164 dense sites moved outward.
   Re-check: the 1-mm sampling guarantees 0.29 mm in general, not 1 mm, and the nearest-triangle
   search was not certified. An exact cube-to-mesh distance (certified candidate triangles and an
   intersection test; `opm.exact_cell_clearance`, tested) puts every final adult cell at least
   1.001 mm out and every cell of the six G3B anatomies at least 0.96 mm out (7 of 1,502 cells at
   0.96-1.00 mm); the rule and the arrays are unchanged and the documentation states this.
2. [x] M2 practical detection scored an injected event by the statistic's maximum in the truth
   window, while false events were the detector's emitted events (local maxima >= refractory
   apart). Injected events are now scored with the emitted events (`detection.event_height`), in
   detection, localization and the ROC.
3. [x] M3 patch-support recovery counted grid sources within the patch radius of its centre, not
   the patch's members: now the members on the inverse grid.
4. [x] M4 censored resamples (a system not reaching 50 % detection, a motion draw not reaching
   the loss level within the tested range) were dropped from the intervals: now kept as censored
   values, with open interval ends where they fall. Re-check: the S50-ratio intervals still dropped
   resamples where neither system reached 50 %, and the 0 / inf sentinels claimed more than is
   known. An S50 outside the tested strengths is now a bound (at or below the weakest, beyond the
   strongest), each resample's ratio an interval, and the 95 % interval takes the outer percentiles
   of the bounds over every resample (`detection.censored_interval`, tested); 37 of the 288 paired
   ratio intervals change (31 visible at two decimals; 29 at 45-70 mm), all in bands with censored
   resamples, re-summarised from the stored simulations at d667a9d.
5. [x] M5 the report builder deleted any existing output directory: now only one it made (marker
   file) or an empty one.
6. [x] m1 Neuromag "magnetometers"/"gradiometers" after the external projection are subsets of the
   jointly projected 306 channels: stated as such.
7. [x] m2 the room-noise synthesis reversed the cross-spectral phase (scipy's CSD convention):
   fixed, with a delayed-channel phase test.
8. [x] m3 the held-out ROC false-event axis was off by one event: the count above each threshold.
9. [x] m4 G1C patch-area scaling: 8.3 dB for the simulated patches (was 8.5).
10. [x] m5 lead-field fingerprints also cover the coil definitions, the source positions and
    normals and the MNE version. Re-check: MNE's standard coil definitions file is now hashed too,
    in the forward-cache key and the fingerprint (tested by moving one integration point). The seven
    lead fields were recomputed without the forward cache: bit-identical to those of 19a8fd2. Their
    sidecars first recorded the commit at the time of writing (a local import); the lead-field
    module now fixes it at import, and the sidecars record bf686e6.
11. [x] m6 report downloads: the nested Fig. 3 files and the motion CSV carry their generating commit.
12. [x] B1 (a public-release blocker, not a defect of the results): the history holds identifying
    metadata and earlier PDFs with licensed fonts; a public release needs a fresh export or a
    rewritten history (owner decision, `docs/release_checklist.md`). Nothing was rewritten.
13. [x] N1 wording: G3B is a geometry/fit experiment and G4 bounded simulated-source recovery, without
    diagnostic, epileptogenic-zone or surgical claims.

Also found while fixing: the first M2 fix reused a variable that the localization injection needs
(caught before any rerun; the scoring now lives in one tested function). In the rerun the adult
localization was first started before the full-resolution lead fields were rebuilt and stopped at
their fingerprint check (nothing written); it ran again after them.

Effect on the results (every OPM-dependent result recomputed at 19a8fd2, `adult-baseline-v3`; G1A,
G3A and the exact-sphere BEM check do not involve the arrays and are kept):
- G2: dense OPM 1.11x [1.08-1.14] Neuromag combined (v2 1.13x [1.10-1.16]); after the 8-term
  projection 1.05x [0.99-1.10], a tie within its interval (v2 1.07x [1.01-1.13]); matched 1.00x.
  Refining the BEM head surface now changes the headline by -0.1 % (v2 -1.6 %), and the cell and a
  point sensor differ by at most 0.5 % (v2 2.8 %): the v3 cells sit farther out; the convergence of the
  refined surface itself remains an assumption.
- G3B: D_adult +0.85 dB (v2 +0.99); every Delta 0.12-0.20 dB larger (+0.59 to +1.25 dB); the
  counterfactual helmet about the laterally centred head -0.21 to +0.01 dB (v2 -0.37 to -0.12).
- G4 adult: dense vs Neuromag at 10-20 mm 8/2 locations, p = 0.04, ratio 1.29 [1.03-1.51] (v2 16/0,
  1.51, with the old scoring and other noise draws; oracle 14/2, p = 0.001); the adult dSPM gain of
  earlier runs did not reproduce (0 mm). Pediatric: superficial ratios 1.33-1.61 (v2 1.33-1.64);
  four localization comparisons survive a correction within anatomy and array (v2 seven, mostly
  others); of 192 localization comparisons 30 have p < 0.05, 29 favouring an OPM array (half of the
  30 for weak, mostly undetected 80-nAm sources).
- Motion: the 1- and 3-dB loss thresholds within a few per cent of v2; the adult's D = 0 thresholds up
  to 32 % lower (its static D after the projections +0.32 to +0.65 dB, v2 +0.51 to +0.80); with the
  artefact modelled the loss stays below 0.17 dB (v2 0.15 dB).

## Equal OPM standoff (2026-10-02/03): `adult-baseline-v4`

Every anatomy's BEM head surface has its vertices on its MRI scalp (A-BEM-CONFORM; decisions log).
Effect on the results (every head-model- and array-dependent result recomputed: the full-resolution
lead fields, G1A, G3A, adult and pediatric detection, the near-mesh and head-surface studies and
motion at ed852b2; G1B, G1C, G2, G3B and every localization at e53bea8; the decomposition at 236068d;
the G2 report, the G3B summaries and figures and the G4 cross-reading summaries at 9512f1b; the
exact-sphere BEM check involves no head model and is kept):
- G2: dense OPM 1.14x [1.12-1.17] Neuromag combined (v3 1.11x); after the 8-term projection 1.12x
  [1.08-1.15] (v3 1.05x, a tie); matched 1.01x, 0.95x projected (v3 1.00x, 0.89x). The shallow ratios
  did not change and the deep ones rose, mostly through the sensitive axes
  (`results/g2/head_surface_effect.json`).
- G3B: D_adult +1.00 dB (v3 +0.85); every Delta 0.11-0.27 dB smaller (+0.44 to +1.07 dB); the
  counterfactual helmet about the laterally centred head -0.38 to -0.16 dB (v3 -0.21 to +0.01).
- G4 adult: dense vs Neuromag at 10-20 mm 11/1 locations, p = 0.004, ratio 1.36 [1.10-1.58] (v3 8/2,
  1.29), and also at 20-30 and 45-70 mm (uncorrected; of the earlier runs only an early v2 run, at 20-30 mm); the matched array's
  deficit at 45-70 mm is not seen (3/6); dSPM errors for 320-nAm patches 5.6 and 2.4 mm smaller
  (p = 0.054 and 0.067; v3 0 mm). Pediatric: the templates' arrays did not change, so their outcomes
  are essentially those of v3; superficial ratios 1.33-1.65; 52 of the 240 localization comparisons have
  p < 0.05, 51 favouring an OPM array.
- Motion: the adult's 1- and 3-dB thresholds move by -33 % to +11 % and its D = 0 thresholds by -25 %
  to +78 % (1 dB in a uniform field without correction 0.017 deg, v3 0.022; its static D after the projections +0.67 to +0.87 dB, v3 +0.32 to +0.65), the templates'
  not at all; with the artefact modelled the loss stays below 0.15 dB.

## Open questions and risks

1. School-aged native anatomy: obtained 2026-10-03 (three children of OpenNeuro ds005234, decisions log), with a
   modelled skull because the dataset's watershed inner skulls lie just below the scalp, and transferred
   fiducials because none come with the data (the talairach.xfm files, shifted by one subject in the
   snapshot like the BEMs, were checked but not used); three individuals of one dataset are not a population, so variability between children is
   shown, not estimated. Their derived figures and tables need the dataset's licence checked (CC0 in its
   metadata, CC BY in its acknowledgement text) before any public release.
6. Resolved in v4 (2026-10-02/03): the adult's BEM head surface lay about 1 mm outside its MRI scalp, the
   templates' on theirs, so the whole-cell clearance rule moved the adult's OPM sensors outward and its stored
   normals tilted their axes; every anatomy's head surface is now on its MRI scalp (A-BEM-CONFORM), every
   array sits at a median 7.0 mm (decisions log, `results/g2/head_surface_effect.json`).
7. Done in v4 from the goal review's open list: the G1C variant without the medial wall and with untruncated
   patches; the bounded localization with `mne.minimum_norm` (next to the study's dSPM), with Neuromag
   magnetometers or gradiometers alone, and with a held-out check of the detector thresholds; a
   measured-spectrum SQUID noise scenario; paired plug-in covariances; a triaxial channel-count control;
   20-mm and fixed-density patches and a yaw placement variant in G3B, with citations for the placement
   (A-G3-PLACE); status lines in the result CSVs and on the G1A replica's summary; a minimum-count rule for
   the G1B patch bins and the segment-only Hilbert variant; `RESUME=1 scripts/run_all.sh`. Still open: a CI
   workflow (owner decision: minutes on a private repository).
2. The TRIUX specification image was not found; values from GOAL.md are used and flagged.
3. OPM device noise: no single verified device specification; a declared 7-30 fT/sqrt(Hz)
   sweep is used instead.
4. Clean-environment smoke test: done on 2026-09-30 with a fresh venv at a commit before 81168f3 (87
   tests then); since then verification used fresh clones with the development venv (133 tests at
   9771f68, 132 pass, 1 skip; 148 tests at db3f0dc, 147 pass, 1 skip; 160 tests at 8bdeb8b, 159 pass, 1 skip). No routine CI workflow exists (owner decision: minutes on a private
   repository); `README.md` gives a smoke command.
5. MEG 2443 is bad in the sample recording (baseline RMS 23x the gradiometer median): exclude it from
   every measured-noise computation (G2 brain-noise calibration, empty-room fit).
