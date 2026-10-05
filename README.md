# OPM vs SQUID MEG: adult-to-pediatric simulation study

Simulation study comparing on-scalp optically pumped magnetometers (OPM) with the Neuromag
(Vectorview/TRIUX-type) SQUID system: an analytical benchmark, adaptations of two published
adult studies, a realistic adult comparison, a pediatric fixed-helmet versus head-adaptive
extension, and epilepsy (interictal spike) detection and localization examples in both.
The milestones, plan, status and decisions log are in `PLAN.md`.

This is a proposed study in development. It does not claim that any author of the reproduced
papers has reviewed or approved it. The OPM advantage is tested, not assumed.

## Status

| Milestone | Label | Status |
|---|---|---|
| G0 audit, provenance, plan | - | done |
| G1A Jas et al. 2026 analytical benchmark | REPRO | done, checked |
| G1B Hunold et al. 2016 depth-orientation spike SNR (MEG part) | ADAPT (+ NEW OPM column) | done, checked |
| G1C Goldenholz et al. 2009 cortical SNR maps (MEG part) | ADAPT (+ NEW OPM extension) | done, checked |
| G2 realistic adult OPM vs Neuromag | NEW | done, checked; frozen as `adult-baseline-v1`, corrected in `adult-baseline-v2`, `adult-baseline-v3` and `adult-baseline-v4` (equal OPM standoff on every head) |
| G3 pediatric extension | NEW (size benchmark REPRO) | G3A size benchmark done; G3B done (24-, 18- and 12-month infant templates, three school-aged children of OpenNeuro ds005234 with a modelled skull, and two scaled-adult size controls), after the adult freeze `adult-baseline-v2`, recomputed with the v3 arrays (`adult-baseline-v3`) and with every head surface on its MRI scalp (`adult-baseline-v4`); the 2-year-template pass checked (fixes re-verified); the 18- and 12-month templates added on 2026-10-01, the school-aged children on 2026-10-03 (v4 only) |
| G4 epilepsy detection and localization | NEW | adult done; pediatric done (three templates, both size controls and three school-aged children); head motion and OPM slippage: bounded secondary extension done; adult and 2-year-template parts checked |
| G5 software, reproduction, report | - | 197 unit tests pass; the report is built from the stored results and published as a GitHub Pages site |

"Checked" means separate checking passes within the project by reviewer agents, with fixes
re-verified (on 2026-10-01 also complete passes by two other models); it is not external peer review.

Adult findings so far (`adult-baseline-v4`), conditional on one adult head (MNE sample subject),
an assumed OPM noise of 15 fT/sqrt(Hz), OPM sensors at the 7-mm standoff with no extra
helmet-to-scalp gap (median sensing-centre height 7.0 mm above the MRI scalp, as on every head
of the pediatric extension) and Neuromag at its measured (not best) head position; details and
caveats in `docs/methods.md` and `results/g2/G2_report.md`:
- With modelled brain noise (calibrated on gradiometers; it predicts 0.73x the measured
  magnetometer brain noise), a dense on-scalp OPM array (208 single-axis sensors, against
  Neuromag's 306 channels) has 1.14x [1.12-1.17] the known-topography detectability of the
  Neuromag system: about 1.6x for sources 10-15 mm below the scalp, falling to 1.05-1.07x below
  35 mm. After the 8-term external-field projection, the study's second headline condition, the
  ratio is 1.12x [1.08-1.15] (higher for 84 % of targets; the projection costs the dense array 3 %
  and Neuromag 0.5 % of their detectability). An OPM array at Neuromag's own 98 sites ties it
  (1.01x), and falls behind after the projection (0.95x [0.90-0.98]). With Neuromag's channel
  count (a triaxial OPM at those sites, 294 channels) it is 1.12x, or 1.05x if the tangential axes
  are twice as noisy. With covariances estimated from 10 or 60 s of data instead of the oracle, the
  dense ratio is 1.23x or 1.16x (Neuromag's 306 channels lose more to estimation). The ratio
  hardly depends on the brain-background level; it is set by signal outside the brain-noise
  subspace.
- The advantage disappears with worse assumptions: 1.01x at an OPM noise of 30 fT/sqrt(Hz),
  0.97x with 30 fT/sqrt(Hz) and a 3-mm scalp gap, 0.92x with a 6-mm gap. Without brain noise,
  Neuromag wins (0.74x). The 3- and 1-layer head models agree (1.14x and 1.13x). The advantage also
  depends on the metric and the depth: by the best single channel (peak-channel SNR) Neuromag is
  ahead (dense 0.90x, Neuromag higher for 72 % of targets; mean-power SNR 1.04x), and after the
  external-field projection the dense array is no longer ahead below 45 mm (0.92-1.00x; the matched
  array is behind from 25 mm, 0.75-0.99x). With Neuromag noise taken from the measured empty-room spectrum
  instead of the brochure values the dense ratio is 1.17x. Not modelled: non-cortical physiological
  fields (cardiac, ocular), which the OPMs would see as magnetometers; the 1/f rise of real OPM
  noise at low frequencies (white noise assumed); and the measured brain noise's spatial pattern,
  which the cortical background matches only in its median gradiometer level (per channel the
  model/measured ratio spans 0.16-1.90).
- Equal standoff (v4): the sample subject's stored outer skin lay about 1 mm outside its MRI scalp,
  so the whole-cell clearance moved most of the adult's OPM sensors outward (v3), and its stored
  normals tilted their axes; with the head surface on the scalp, as the infant templates' already
  are, the shallow ratios are unchanged (1.6x at 10-15 mm) while deep sources gained (projected
  below 45 mm 0.92-1.00x, v3 0.69-0.76x): the dense array now also covers the lower occipital scalp
  that v3's outer skin excluded, and the projection costs it less. A decomposition
  (`results/g2/head_surface_effect.json`) attributes the projected-condition change (1.05x -> 1.12x)
  mainly to the sensor axes (the stored normals were nearly radial: 1.07x -> 1.12x) and partly to the
  sites (standoff and coverage: 1.05x -> 1.07x; the 12 lower occipital sites alone about 2 %); the
  conductor surface itself changes nothing.
- Simulated interictal spikes: at 10-20 mm depth the dense array needs about a quarter less source
  strength for 50 % detection (34 [28-51] vs 47 [36-63] nAm; the intervals overlap, the paired
  strength ratio is 1.36 [1.10-1.58]); 11 of the 18 locations favour the OPM, 1 Neuromag, 6 tie
  (p = 0.004, uncorrected); with an oracle detector that knows the source and onset, 15 favour the
  OPM and none Neuromag (p = 0.00006). The practical detector knows the three simulated spike
  morphologies, so absolute sensitivities and false-event rates are optimistic; the paired
  comparison is less affected. The superficial advantage holds in every run, its size varies (v3
  8/2, ratio 1.29; v2 16/0, 1.51). In v4 the dense array is also favoured at 20-30 mm (9/1,
  p = 0.014) and 45-70 mm (9/0, p = 0.004); the equal-standoff change raised its G2 detectability
  at the G4 locations in every band, most in the deepest, but 30-45 mm shows no detection difference.
  With 16 such tests per run, other noise draws (the channel count changes the
  random stream) and earlier runs without a deeper difference (v3: p = 0.19 and 0.50), these are
  reported, not established. The matched array is not ahead in any band (10-20 mm: 7/4, p = 0.40)
  and no longer behind at 45-70 mm (3/6, p = 0.40; v3 1/9, p = 0.016). v1's claim that the dense
  array was favoured in every depth band (14/0, 10/2, 12/2, 13/0) is superseded. The 72 locations
  are frontal-heavy (24 frontal, 2 occipital).
- Bounded localization (24 locations, one event each): dipole errors are similar across arrays
  (about 4-6 mm for strong focal spikes, limited by a 2-mm/2-deg coregistration error; no dipole fit
  failed, and dipoles more than 30 mm off are mostly undetected weak events (adult: 81 of 111),
  `results/g4/G4_fit_failures_report.md`). For extended 320-nAm sources both OPM arrays are detected
  and dSPM-localized within 10 mm at more locations than Neuromag (8/0 and 7/0, p = 0.008 and 0.016,
  uncorrected; not surviving a correction over the 20 comparisons per array), with dSPM errors 5.6
  and 2.4 mm smaller (p = 0.054 and 0.067); computed with MNE's own `mne.minimum_norm` the
  differences are 1.0 and 2.3 mm (p = 0.22 and 0.073). Earlier runs gave 0 mm (v3) and about 4-5 mm
  (v2): not robust. Differences not detected are not excluded.
- `adult-baseline-v1` reported 1.21x; that value was inflated by OPM cells reaching into the
  coarse BEM head surface, a numerical error corrected in v2 (see PLAN.md). v2 (1.13x) kept the
  cells' integration points outside the head surface but not their corners; v3 (after later checks)
  kept the whole cell outside (1.11x), which with the sample's outer skin about 1 mm outside its
  scalp moved most sensors about 1 mm out; v4 puts the outer skin on the scalp (equal standoff
  with the children), which leaves 208 dense and 98 matched sites (1.14x).

Pediatric findings (G3B, NEW), conditional on one adult head, three average infant templates of
one database (24, 18 and 12 months; O'Reilly et al. 2021), three individual school-aged children of
one dataset (7.8-8.7 years; OpenNeuro ds005234, Fadeev et al. 2024; their skull modelled because the
dataset's segmentation failed, their fiducials transferred from the adult) and two scaled copies of
the adult; the same Neuromag helmet, sensors and noise for every head; OPM arrays refitted to each head with the
adult rules; details in `docs/methods.md` section 10 and `results/g3b/G3B_report.md`:
- With the child raised to 20-mm contact with the top of the fixed helmet, the dense OPM array's
  known-topography detectability relative to Neuromag combined (D) rises from +1.00 dB in the adult
  to +1.50 dB (school-age-size control), +2.13 dB (2-year-size control), +1.84 dB (2-year
  template), +1.89 dB (18-month template) and +2.04 dB (12-month template): Delta = D_child -
  D_adult = +0.44 [+0.34, +0.56], +1.07 [+0.91, +1.18], +0.73 [+0.46, +0.98], +0.88 [+0.54, +1.05]
  and +0.96 [+0.61, +1.51] dB (for the scaled controls the vertex-wise difference; for the templates
  the adult-area-weighted median over parcels of the parcel difference, which is not the difference of
  the two medians, e.g. +0.73 vs +0.85 dB for the 2-year template). Every head's OPM sensors sit at
  the same median 7.0 mm above its MRI scalp (v4; in v3 the adult's sat 0.8 mm farther out and
  Delta was 0.11-0.27 dB larger). Against the gradiometers or magnetometers alone, after the
  external-field projection (of all 306 channels jointly), for the matched-site OPM array and for
  extended sources (5-, 10- and 20-mm patches) the sign is the same (for the school-aged children
  below, except child C after the projection).
- In the templates and size controls the gain comes from the fixed helmet's fit: left at the adult's
  ear-line position Delta is +1.14 to +2.50 dB; laterally centred or at true 18-mm contact it stays
  +0.44 to +1.09 dB; in a counterfactual helmet scaled with the head, centred laterally, it is -0.21,
  -0.38, -0.30, -0.16 and -0.18 dB (without the lateral centring the off-centre 24- and 18-month
  templates kept +0.24 and +0.18 dB). With the background fixed per unit cortical area, both systems'
  detectability rises in these smaller heads and the on-scalp OPM's rises more (vertex-wise +1.10 and
  +1.67 dB vs Neuromag +0.65 and +0.54 dB in the scaled controls); with a helmet that fits (scaled
  with the head, laterally centred), the SQUID gains at least as much. The children's arrays also
  have fewer OPM sites (174, 155, 151, 157 and 144; school-aged children 155, 153 and 167; adult
  208); at an equal site count Delta would be larger.
- The child-minus-adult difference of the median D stays positive for OPM noise 7-30 fT/sqrt(Hz),
  background variance x0.5 or x2 and a 1-layer head model (for the school-aged children except at 7
  fT/sqrt(Hz), -0.15 to -0.03 dB); at 30 fT/sqrt(Hz) the adult's D is -0.05 dB (Neuromag slightly
  ahead), the templates' +0.61, +0.61 and +0.77 dB and the school-aged children's +0.55, +0.72 and
  +0.27 dB. In the scaled controls the gain is large near the
  surface (vertex-wise +2.40 dB at 10-15 mm for the 2-year size) and smallest at 40-60 mm. In the
  templates the within-stratum gain is smallest at 20-50 mm, larger near the surface and deeper (most
  deep in the 24- and 18-month templates, at 10-15 mm in the 12-month one) and larger for radial
  sources, but much of their pooled gain reflects their shallower cortex (reweighted to the adult's depth mix, +0.32, +0.53 and +0.29
  instead of +0.85, +0.90 and +1.05 dB at target level); it is also regionally asymmetric, with the
  templates off-centre in the helmet. Top contact is the adult's second-lowest of 12 source-blind
  placements (mid-family for the templates and size controls; 4th, 6th and 3rd from the lowest for
  the school-aged children): against each head's family median the child-minus-adult differences of
  the medians are 0.14-0.20 dB smaller (school-aged children 0.05-0.17 dB).
- Native school-aged heads gain less than the infant templates and, at a similar head circumference,
  less than the size-only control (children A and C; their intervals overlap or touch the control's): Delta +0.30 [+0.11, +0.38], +0.41 [+0.28, +0.54] and +0.17 [-0.07,
  +0.34] dB for children A, B and C (head circumference 520, 486 and 535 mm; the school-age size
  control, 525 mm, +0.44 [+0.34, +0.56] dB; child B's smaller head gains about as much), after the
  external-field projection +0.03, +0.26 and -0.19 dB. Their cortex is nearly adult-sized
  (1,666-1,869 cm^2 of usable cortex, adult 1,878, school-age size control 1,470), so with the
  background fixed per unit area Neuromag's brain noise stays near the adult's while the closer OPM
  sees more of it; reweighted to the adult's depth mix their gain is -0.05 to +0.05 dB, and within
  depth strata it is near zero or negative down to 30 mm and positive only deeper. A scaled adult
  therefore likely overstates the gain of school-aged heads of its size (point estimates); the 2-year size control and the
  infant templates, whose cortices are smaller, may too. In a helmet scaled with the head, about the
  laterally centred head, Delta is -0.76, -0.59 and +0.04 dB (centred -0.74, -0.23 and +0.56 dB;
  child C's helmet barely shrinks).
- At 100 nAm (detectability >= 5, an operational threshold) both systems reach 66 % of the
  adult's, 63-65 % of the school-aged children's and 75-79 % of the templates' usable cortex, the OPM
  alone a further 2 %, 4 % and 5 %, and the SQUID alone at most 0.14 % (0.75-1.18 % in the
  school-aged children). At a fixed current density of 0.5 nAm/mm^2, 10-mm patches (about
  140 nAm) reach it at 43 % (Neuromag) and 54 % (dense OPM) of the adult's patch centres and at
  66-77 % and 76-85 % on the templates.
- Simulated spikes (same detectors and seeds as the adult): at 10-20 mm more locations favour the
  full OPM array in the adult and every smaller head (strength for 50 % detection
  at 10-20 mm, Neuromag combined / dense OPM: adult 47/34; templates and size controls 53/33, 61/37,
  43/30, 36/27 and 43/29 nAm, paired strength ratio 1.33-1.65 against the adult's 1.36, 13-17 of 18
  locations favour the OPM in each; the three school-aged children 61/48, 65/55 and 52/36 nAm, ratio
  1.18-1.43, 11-14 locations, p = 0.004-0.022, none of the three surviving the correction over the 24
  detection comparisons); deeper, location-level differences appear for the adult (20-30 and 45-70
  mm), the 2-year size control (30-70 mm), the 12-month template (30-45 mm) and the 18-month template
  (45-70 mm), not for the school-aged children, uncorrected and none surviving a correction, so G3B's
  deeper gain is mostly not resolved at this sample size. The frozen thresholds give unequal held-out
  false-event rates (0.4-1.55 per minute); with every detector set to 1 per minute on the held-out
  null (`results/g4/G4_matched_rate_report.md`) the superficial result is unchanged (ratios
  1.21-1.63), of the children's deeper ones the 12-month template's 30-45 mm remains (p = 0.018,
  uncorrected) and one school-aged child's 20-30 mm appears (p = 0.032, uncorrected). Localization:
  dipole errors similar; dSPM of strong focal spikes 3.6-7.0 mm better with the dense OPM in both size
  controls and the 2-year and 18-month templates (p = 0.0002 to 0.035, uncorrected; MNE's own dSPM
  agrees in direction, -5.0 to 0.0 mm), not in the 12-month template or the adult (0 mm), and 2.9-4.7
  mm in the school-aged children without a clear difference (p = 0.071-0.25; MNE's dSPM -4.9 to -5.5
  mm, p = 0.009-0.020). Eight pediatric localization comparisons survive a Bonferroni correction
  within their anatomy and array (20 each), none of them in the school-aged children, and none across
  all eight children (320 comparisons; across the five earlier children two did: the school-age
  control's dense dSPM of strong focal spikes and the 12-month matched array's dSPM of weak, mostly
  undetected 80-nAm patches); which ones survive varies between runs, the direction does not: of the
  360 localization comparisons in the nine anatomies, 65 have p < 0.05 (uncorrected; about 18 would
  be expected by chance if they were independent, which they are not), 64 of them favouring an OPM
  array; 28 concern weak, mostly undetected 80-nAm sources.

Scope: G3B is a geometry and helmet-fit experiment on average templates, three individual
school-aged children and scaled copies of one adult, and the G4 examples are bounded recovery of
simulated sources. Neither shows a diagnostic benefit, identifies an epileptogenic zone or
establishes a surgical benefit.

Head motion and OPM slippage (G4, NEW, a bounded secondary extension; `docs/methods.md` section 12,
`results/g4/G4_motion_report.md`), adult and the 24- and 12-month templates, declared conditions:
- A sustained head displacement in the fixed helmet costs Neuromag, analysed with the template of
  the reference position, 0.3-0.6 dB at 5 mm and 1.4-1.7 dB at 10 mm (median detectability), and up
  to 2.3 dB for a 10-deg pitch or roll of a template; with the displaced geometry known (ideal
  movement compensation) at most 0.4 dB. The head-mounted array is unaffected by head displacement;
  if its cap slips by 3 deg (about 4-6 mm at the sensors) it loses 0.1-0.5 dB uncompensated and at
  most 0.1 dB with the slip known.
- In-band head rotation moves the OPM array through the room's residual static field (the SQUIDs
  are fixed and have no such term). If the artefact is not modelled, the dense array loses 1 dB at
  about 0.02 deg RMS of in-band rotation in a 1-nT residual field without correction, and at 0.4-0.5 deg
  after a homogeneous or 8-term field projection with 1-deg/1-% sensor calibration errors (0.13-0.15
  deg with 3 deg/3 %); the thresholds scale inversely with the field (0.04-0.05 deg in 10 nT), and a
  homogeneous projection alone leaves the gradient term (1 dB at 0.14-0.21 deg in a 1-nT/m
  gradient). If the artefact is modelled (oracle covariance), the loss stays below 0.15 dB.
  These are bounds on two mechanisms, not a motion-robustness result; real head-motion statistics,
  sensor dynamic range and gain changes are not modelled.

Labels: REPRO = reproduction with the paper's definitions; ADAPT = adaptation where original
data or details are unavailable; NEW = new experiment or study choice.

## Layout

| Path | Content |
|---|---|
| `src/opmsquid/` | package: sphere model, sensors (Neuromag, OPM), anatomy (adult, infant templates, school-aged children with a modelled skull, scaled controls), FreeSurfer readers without nibabel (`fsio`), forward models (cached), noise, metrics, background, environment, paper-specific modules (`hunold`, `goldenholz`), G2 comparison, pediatric helmet placements (`pediatric`) |
| `scripts/` | one driver per milestone (`g1a_*`, `g1b_*`, `g1c_*`, `g2_*`, `g3a_*`, `g3b_*`, `g4_*`), full-resolution forward jobs, Fig. 6 digitiser, BEM accuracy studies, the school-aged children's fetch, preparation and checks, `run_all.sh` |
| `configs/` | paper and study configurations, every value tagged as printed, chosen or digitised |
| `tests/` | `unittest` suite (physics, units, metrics, geometry, noise, statistics, report builder) |
| `site/` | templates and style of the local report (`scripts/build_site.py`; output in `site/_build/`, not committed) |
| `results/` | figures, CSV and JSON per milestone (JSON records the code commit and versions) |
| `docs/` | methods, provenance register, audit, structured literature extractions |
| `legacy/` | earlier Fig. 3 replication work (science unchanged; open fonts and relative paths since 2026-10-01) |

## Reproduce

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -c "import mne; mne.datasets.sample.data_path(path='data/external')"  # ~1.6 GB download
.venv/bin/python -c "import mne; [mne.datasets.fetch_infant_template(a, subjects_dir='data/external/infant_subjects') for a in ('2yr', '18mo', '12mo')]"  # ~1.15 GB (G3B, pediatric G4)
.venv/bin/python scripts/fetch_school_subjects.py  # ~125 MB: three school-aged children of OpenNeuro ds005234 (G3B, pediatric G4)
bash scripts/run_all.sh     # tests, then every milestone in order (~11 h on a 10-core laptop: adult ~2 h, eight pediatric G4 runs ~1 h each)
```

Smoke run (about 15 min: the unit tests, which skip the full-resolution lead-field check without
the local cache, then the analytical G1A benchmark, which needs no data):

```bash
.venv/bin/python -m unittest discover -s tests -t . && .venv/bin/python scripts/g1a_jas_benchmark.py
```

Each driver writes `results/<milestone>/` and records the code commit it ran in its JSON.
Commit hashes in result files identify the code version in this repository; `docs/commit_map.tsv` maps the hashes
recorded before 2026-10-04 to the current history.

The sample data must end up in `data/external/MNE-sample-data` (`src/opmsquid/paths.py`;
override with `OPMSQUID_DATA`). Caches go to `cache/` (`OPMSQUID_CACHE`).

## Local report

`.venv/bin/python scripts/build_site.py` builds a static report into `site/_build/` (git-ignored;
open `site/_build/index.html`) from the committed result files, in the order the goal asks for:
adult benchmarks, realistic adult results, pediatric extension, epilepsy, then methods,
parameters and reproduction. Every number is read from `results/`; the downloads list each file's
size, SHA-256 and the code commit recorded in it; the build fails on any broken link. The published report,
https://khan-laiba.github.io/opm-squid-meg/, is this build, committed to the `gh-pages` branch by
`scripts/deploy_pages.py` and checked with `scripts/check_live_site.py`.

## Data and privacy

- The infant templates (24, 18 and 12 months; O'Reilly et al. 2021, from the Neurodevelopmental
  MRI Database of Richards et al. 2016; LGPL-2.1 repository) are not committed. Figures and tables
  derived from them (every `results/g3b` figure that shows a template, `results/g3b/g3b_targets_infant*`,
  `results/g4/*infant*` and the motion results) cite both papers. Publishing these derived results
  was decided on 2026-10-04 (the owner's instruction to publish the report with all its analyses;
  `docs/release_checklist.md` item 3): the templates are distributed publicly by their authors
  (J. E. Richards, who created the database, is a co-author of O'Reilly et al. 2021) under LGPL-2.1
  through MNE-Python's `fetch_infant_template`, and only derived results are published.
- The school-aged children (OpenNeuro ds005234, Fadeev et al. 2024; individual, de-identified MRIs of
  typically developing children) are not committed (`data/external/school_subjects/`, fetched by
  `scripts/fetch_school_subjects.py` and listed with S3 object versions, sizes and SHA-256 in
  `configs/school_subjects_manifest.json`). Results derived from them (`results/g3b/g3b_targets_child*.csv`,
  `school_subjects_preparation.json`, `school_anatomy_checks.json`, the G3B summary, report and figures,
  `results/g4/*child*` and the pediatric G4 cross-reading summaries and reports) are derived quantities
  only; the dataset's metadata say CC0 and its acknowledgement text CC BY: the report and the supplementary
  pages cite the dataset and its paper.
- Not committed: the reference PDFs, the MNE sample data and anatomy (`data/`), computed caches
  (`cache/`). Results contain only derived quantities of the public datasets.
- Release preparation (`docs/release_checklist.md`): every figure is set in open fonts (DejaVu
  Sans, STIX); the published Jas et al. Fig. 3 raster used by the legacy verification is CC-BY 4.0
  and carries its attribution; output summaries no longer record local absolute paths; every
  commit and tag records Laiba Khan with GitHub's no-reply address, and no file in the history
  holds an e-mail address. Licences (2026-10-04): the code is under the MIT licence (`LICENSE`); the
  report text, documentation, figures and result files under CC BY 4.0 (`LICENSE-CONTENT.md`).
  Still open: the earlier font-embedding PDFs that remain in the history, the redistribution of fsaverage-derived
  figures, and two documents that name local folders (the template-derived results: decided
  2026-10-04, above).
- Commit hashes: on 2026-10-04 the history was rewritten to correct the recorded author identity
  (file contents and dates did not change, every commit hash did): commit hashes recorded before
  that date, in result files and documents, refer to the original history; `docs/commit_map.tsv`
  maps each of them to the current commit.
- The repository is public (since 2026-10-04) and the report is published at
  https://khan-laiba.github.io/opm-squid-meg/.
