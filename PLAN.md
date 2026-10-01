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
| G1A Jas analytical benchmark | done; independently reviewed (approve with notes; notes addressed) | `scripts/g1a_jas_benchmark.py` -> `results/g1a/`; Eq. 1 vs Sarvas 2-D maximum 4.7e-15, vs MNE sphere 5.4e-8; d_eq(eta 3) = 27.665 mm |
| G1B Hunold depth-orientation spikes | done; independently reviewed (approve with notes; notes addressed); rerun with the 4-mm source rule and a realization-averaged Fig. 6 calibration | `scripts/g1b_hunold.py` -> `results/g1b/`; calibrated background (scalar 0.47; one realization alone 0.44-0.52), p2p: bin means 0.78-0.91x the paper (0.68-1.01x over the calibration range), noisy p2p 0.97-1.04x, r 0.91-0.97, 2.5-classification agreement 80-93 %, GM-MM sign agreement 100 % |
| G1C Goldenholz cortical SNR maps | done; independently reviewed (approve with notes; notes addressed); rerun with the 4-mm source rule | `scripts/g1c_goldenholz.py` -> `results/g1c/`; s_s 1.76 nAm at 4,000 sources (paper 1.6-1.9); focal median -22.1 dB, 56 % inside the paper's -29/-19 dB range; deep medial cortex darkest |
| G2 realistic adult OPM-Neuromag comparison | done; independently reviewed; `adult-baseline-v1` (1.21x) corrected in `adult-baseline-v2` (near-surface BEM error, see below) | `scripts/g2_adult_comparison.py` -> `results/g2/` (`G2_report.md`); with modelled brain noise the dense 212-site OPM array is 1.13x [1.10-1.16] Neuromag combined (1.01-1.28x for OPM noise 30-7 fT/sqrt(Hz), 0.93-1.13x jointly with a 0-6 mm scalp gap; 3- and 1-layer BEM agree), the matched 97-site array 1.00x (a tie); 1.4-1.7x for sources within 20 mm of the scalp, 1.02-1.05x below 35 mm; without brain noise Neuromag wins (dense 0.76x); all adult results rerun at 81168f3 |
| G3 pediatric extension | G3A done (REPRO size benchmark); G3B done after the adult freeze (24-, 18- and 12-month infant templates, adult scaled to school-age and 2-year head size; fixed Neuromag helmet at source-blind placements vs refitted OPM; counterfactual helmet) | `scripts/g3a_jas_size_benchmark.py` -> `results/g3a/`; `scripts/g3b_pediatric_helmet.py` -> `results/g3b/` (`G3B_report.md`; computed at cb1b8a9, redrawn at fd40dfe): dense OPM vs Neuromag combined at top contact, D_adult +0.99 dB, Delta +0.47 (school-age size), +1.13 (2-year size), +0.73 (2-year template), +0.85 (18 months), +1.03 dB (12 months); with a helmet scaled with the head about the laterally centred head -0.17, -0.30, -0.37, -0.29, -0.12 dB (the gain is mainly the fixed helmet's fit); robust to OPM noise 7-30 fT/sqrt(Hz), background x0.5/x2, 1-layer BEM |
| G4 epilepsy detection and localization | adult (rerun at 81168f3) and pediatric (three templates, both size controls) done; motion and slippage bounded extension done (5bfe1c0) | `scripts/g4_epilepsy_adult.py`, `scripts/g4_localization.py` -> `results/g4/`; 50 % detection at 1 false event/min: dense OPM 35 vs Neuromag combined 53 nAm at 10-20 mm (paired ratio 1.51 [1.25-1.66]; 16/0 locations, p < 0.001), 271 vs 290 nAm at 45-70 mm; no location-level advantage established deeper than 20 mm (the 20-30 mm result varied between runs); matched: 9/3 locations at 10-20 mm (p = 0.04 uncorrected; 0.14 in an earlier run) and 0/7 at 45-70 mm (p = 0.016 uncorrected; same direction earlier); localization: ECD similar (~4-5 mm), dSPM of 320-nAm patches ~5 mm better with OPM (not established: p = 0.14 / 0.046 uncorrected); pediatric (`scripts/g4_epilepsy_pediatric.py`, `G4_pediatric_report.md`): 10-20 mm strength ratio Neuromag/OPM 1.33-1.64 as in the adult (1.51), deeper differences only for the 2-year size control and the 18-month template (uncorrected); motion (`scripts/g4_motion.py`, `G4_motion_report.md`): uncompensated head displacement costs Neuromag 1.4-1.7 dB at 10 mm, a 3-deg OPM cap slip 0.1-0.5 dB; un-modelled in-band rotation costs the OPM 1 dB at ~0.02 deg RMS in a 1-nT field (no correction) or ~0.4 deg (field projection, 1 deg/1 % calibration), thresholds inversely proportional to the field |
| G5 repository, tests, report | private repository `khan-laiba/opm-squid-meg` (no Pages, verified 2026-10-01); 120 unittest tests pass (also in a clean clone of 5bfe1c0); `scripts/run_all.sh`; local report (`scripts/build_site.py` -> `site/_build/`, link-checked, not deployed); clean-environment smoke test passed (2026-09-30, at a commit before 81168f3, 87 tests then); release prepared, not executed | `docs/release_checklist.md`: figures in open fonts, no local paths in output summaries; owner decisions left: licence, git author metadata, template/fsaverage-derived figures, the CC-BY raster, two documents naming local folders, visibility and Pages |

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
| School-aged anatomy (FreeSurfer surfaces, BEM surfaces, scalp) | G3 | **missing** (owner decision 2026-09-30: route (c) + (d)). FreeSurfer is not installed; MNE only packages infant templates (up to 2 years, Neurodevelopmental MRI Database via `mne.datasets.fetch_infant_template`). Routes: (a) age-specific average templates of the Neurodevelopmental MRI Database (registration required); (b) an OpenNeuro school-aged dataset that ships FreeSurfer outputs (e.g. fMRIPrep `sourcedata/freesurfer`); (c) the 2-year infant template as a young-child anatomy (**downloaded**, `data/external/infant_subjects/ANTS2-0Years3T`), and since 2026-10-01 the 18- and 12-month templates of the same series (**downloaded**); (d) scaled adult surfaces as a labelled size-only control (**used**: school-age and 2-year size). A school-aged native anatomy and individual children (between-child variability) remain open. |
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
3. [x] SNR as printed; p2p primary (Fig. 6 digitisation), peak and noisy peak as variants.
4. [x] Absolute level: the paper's Fig. 6 baselines are 0.47x (MM) / 0.42x (GM) our expected
   baselines (20 realizations); both levels reported (U-HU-bglevel). Calibrated, p2p: strong
   bins 0.87-0.92x the paper, weak bins 0.73-0.90x; noisy p2p 0.97-1.04x (the paper lies between
   the two numerators). GM > MM superficially, convergence with depth: reproduced.
5. [x] OPM (matched 97 sites, NEW): brain noise only, OPM > MM for superficial sources (1.4x at
   20-25 mm), equal near 40-45 mm, <= MM deeper (0.9x), < GM everywhere (0.70-0.80x, best
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
   controls. A school-aged native anatomy is still missing, and with average templates of one
   database no claim about anatomical variability is made.
3. [x] Centred, top- and back-contact and bounded translated/rotated placements chosen
   source-blind; regional gaps (helmet selections), source-to-sensor distances, OPM coverage and
   packing, noise composition; counterfactual helmet scaled with the head (labelled mechanistic
   control).
4. [x] D_child, D_adult, Delta for each SQUID comparator, vertex-wise (scaled controls) and per
   parcel and depth/orientation stratum (template), area-weighted, parcel-bootstrap intervals,
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
   fd40dfe): the four anatomies of the reviewed pass reproduce to 1e-9 (188 secondary bootstrap
   bounds moved with the shared random stream); Delta +0.85 [+0.49, +1.27] and +1.03 [+0.66, +1.44]
   dB; the depth checks of the earlier hand computation are now computed for every template
   (`template_depth_checks`).

## G4 — epilepsy (NEW)

1. [x] (adult) IED-like events across regions/depths/orientations/extents and strengths.
2. [x] (adult) Oracle vs practical detector; thresholds calibrated on independent null data and
   frozen; sensitivity vs false events/min; paired comparisons with the location as the unit.
3. [x] (adult) Bounded ECD and MNE/dSPM localization with off-grid truth and registration/model
   mismatch (coregistration draws shared by all arrays); nondetections and errors reported
   separately; GOF descriptive only (whitened GOF depends on the channel count).
4. [x] Pediatric runs (`scripts/g4_epilepsy_pediatric.py`): the adult detection and localization studies
   unchanged on the three templates and both size controls (Neuromag at top contact, refitted OPM
   arrays); comparison in `results/g4/G4_pediatric_report.md` (regenerated at fd40dfe).
5. [x] Motion and slippage, a bounded secondary extension (`scripts/g4_motion.py`,
   `src/opmsquid/motion.py`, `results/g4/G4_motion_report.md`; methods section 12).

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

## Open questions and risks

1. School-aged native anatomy: not available (owner decision 2026-09-30: the 2-year infant template and
   scaled-adult size controls; 2026-10-01: the 18- and 12-month templates added); a claim about anatomical
   variability needs individual children, not averages of one database.
2. The TRIUX specification image was not found; values from GOAL.md are used and flagged.
3. OPM device noise: no single verified device specification; a declared 7-30 fT/sqrt(Hz)
   sweep is used instead.
4. Clean-environment smoke test: done (passed, 2026-09-30).
5. MEG 2443 is bad in the sample recording (baseline RMS 23x the gradiometer median): exclude it from
   every measured-noise computation (G2 brain-noise calibration, empty-room fit).
