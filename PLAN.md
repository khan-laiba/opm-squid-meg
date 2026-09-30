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
| G0 audit, provenance, plan | done except the first push (private repository) | this file; `docs/audit.md`; `docs/literature/`; `docs/provenance_register.md` |
| G1A Jas analytical benchmark | done; independently reviewed (approve with notes; notes addressed) | `scripts/g1a_jas_benchmark.py` -> `results/g1a/`; Eq. 1 vs Sarvas 2-D maximum 4.7e-15, vs MNE sphere 5.4e-8; d_eq(eta 3) = 27.665 mm |
| G1B Hunold depth-orientation spikes | done; independently reviewed (approve with notes; notes addressed, rerun) | `scripts/g1b_hunold.py` -> `results/g1b/`; calibrated background, p2p: bin means 0.84-0.97x the paper (0.69-1.03x over the calibration range), r 0.93-0.97, 2.5-classification agreement 86-95 %, GM-MM sign agreement 100 % |
| G1C Goldenholz cortical SNR maps | done; independently reviewed (approve with notes; notes addressed, rerun) | `scripts/g1c_goldenholz.py` -> `results/g1c/`; s_s 1.78 nAm at 4,000 sources (paper 1.6-1.9); focal median -22.0 dB, 55 % inside the paper's -29/-19 dB range; deep medial cortex darkest |
| G2 realistic adult OPM-Neuromag comparison | done, awaiting independent review, then freeze (`adult-baseline-v1`) | `scripts/g2_adult_comparison.py` -> `results/g2/`; with brain noise the dense OPM array is 1.17x Neuromag combined (1.02-1.34x over the sensitivity analyses), the matched array 1.01x; the advantage is 1.4-1.6x for sources within 20 mm of the scalp |
| G3 pediatric extension | not started; needs pediatric anatomy (see Inputs) | |
| G4 epilepsy detection and localization | not started | |
| G5 repository, tests, report | repository initialised locally | |

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

## Inputs and availability

| Input | Needed for | Status |
|---|---|---|
| Jas et al. 2026 preprint (PDF) | G1A, G3 size benchmark | available locally (not committed) |
| Hunold et al. 2016 (PDF) | G1B | available locally (not committed) |
| Goldenholz et al. 2009 (PDF) | G1C | available locally (not committed) |
| TRIUX specification image | G2 noise and geometry references | **not found on disk**; values transcribed in GOAL.md are used and flagged |
| MRN Neuromag page (T3 coils 3014/3024, shielded-room floor 5-7 fT) | G2 | read 2026-09-30 |
| MNE implementation docs (coil definitions, frames), MNE 1.13.2 | all | read 2026-09-30 |
| Neuromag 3-D geometry, `dev_head_t`, head-MRI trans | G2 | MNE sample data (MGH Vectorview: 204 grads coil 3012 T1, 102 mags coil 3024 T3) |
| Measured SQUID noise | G2 measured-noise scenario | MNE sample `ernoise_raw.fif` (empty room) and baseline covariance |
| Adult anatomy | G1B, G1C, G2 | MNE `sample` subject (3-layer BEM surfaces, oct-6 source space); fsaverage (sibling folder, read-only) |
| School-aged anatomy (FreeSurfer surfaces, BEM surfaces, scalp) | G3 | **missing**. FreeSurfer is not installed; MNE only packages infant templates (up to 2 years, Neurodevelopmental MRI Database via `mne.datasets.fetch_infant_template`). Candidate routes (owner decision + download approval at G3): (a) age-specific average templates of the Neurodevelopmental MRI Database (registration required); (b) an OpenNeuro school-aged dataset that ships FreeSurfer outputs (e.g. fMRIPrep `sourcedata/freesurfer`); (c) the 2-year infant template as an additional young-child anatomy; (d) scaled adult surfaces as a labelled size-only control only. |
| Original Hunold/Goldenholz participant data and recordings | exact reproduction | unavailable; G1B and G1C are adaptations |
| Jas SEF recordings | experimental validation | unavailable; used as context only |

## G0 — audit, provenance, plan

- [x] Goal files, condensed goal, git repository, pinned `requirements.txt`.
- [x] Legacy work moved to `legacy/` with a status table (`legacy/README.md`).
- [x] Full-text extractions: `docs/literature/{jas2026,hunold2016,goldenholz2009}.md` (each independently verified).
- [x] `docs/audit.md`: validated vs unverified legacy claims; reusable components.
- [x] `docs/provenance_register.md`: paper-reported (J, HU, GO), hardware (HW), recovered (R), new assumptions (A), decisions (D); unresolved ambiguities listed.
- [x] `docs/methods.md` (living document).
- [ ] Private repository created, first push, tag `g0`.

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
   stratified to the paper's per-bin counts; 20-mm^2 patches (2895 grow).
2. [x] Background: 10 % of vertices, band-limited, stationary (edge-transient defect fixed),
   +/-10 nAm per dipole; one realization shared by sources and arrays.
3. [x] SNR as printed; p2p primary (Fig. 6 digitisation), peak and noisy peak as variants.
4. [x] Absolute level: the paper's Fig. 6 baselines are 0.43x (MM) / 0.39x (GM) ours; both levels
   reported (U-HU-bglevel). Calibrated: strong bins 0.96-1.04x the paper; weak bins
   0.83-0.92x (paper floor near 1). GM > MM superficially, convergence with depth: reproduced.
5. [x] OPM (matched 99 sites, NEW): brain noise only, OPM > MM for superficial sources, <= MM
   below ~45 mm, < GM everywhere; intrinsic sensor noise (3.5-30 fT/sqrt(Hz)) is negligible
   against this background.
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
6. Freeze the adult baseline as git tag `adult-baseline-v1`.

## G3 — pediatric extension (NEW; size benchmark REPRO)

1. REPRO Jas Table 1 / Fig. 5 size-following benchmark (s = h + 18 mm).
2. Fixed adult Neuromag helmet vs OPM refit per head; physically scaled pediatric anatomy
   (school-aged first) - **input missing**; scaled adult only as a labelled size-only control.
3. Neutral and bounded translated/rotated SQUID placements chosen source-blind; regional
   gaps, coverage, sensitivity; counterfactual reduced-mismatch control (labelled).
4. D_child, D_adult, Delta for each SQUID comparator, per homologous region and strata, with
   uncertainty; show unmatched strata.

## G4 — epilepsy (NEW)

1. IED-like events across regions/depths/orientations/extents and strengths.
2. Oracle vs practical detector; thresholds calibrated on independent null data and frozen;
   sensitivity vs false events/min.
3. Bounded ECD and MNE/dSPM localization with off-grid truth and registration/model mismatch;
   nondetections, failed fits and errors reported separately.

## G5 — software, reproduction, report

Modular package `src/opmsquid`, scripts per milestone, `configs/`, `tests/` (unittest; the
venv has no pytest), cached forward models keyed by input hashes, fixed seeds, resumable runs,
clean-environment smoke test (needs a package download; ask first), local static report in
`site/` (jinja2), never deployed without approval.

## Open questions and risks

1. School-aged anatomy source (G3) - owner decision needed.
2. The TRIUX specification image was not found; values from GOAL.md are used and flagged.
3. OPM device noise: no single verified device specification; a declared 7-30 fT/sqrt(Hz)
   sweep is used instead.
4. Clean-environment smoke test requires installing packages into a fresh venv (download).
5. MEG 2443 is bad in the sample recording (baseline RMS 23x the gradiometer median): exclude it from
   every measured-noise computation (G2 brain-noise calibration, empty-room fit).
