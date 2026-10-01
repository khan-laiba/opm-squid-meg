# OPM vs SQUID MEG: adult-to-pediatric simulation study

Simulation study comparing on-scalp optically pumped magnetometers (OPM) with the Neuromag
(Vectorview/TRIUX-type) SQUID system: an analytical benchmark, adaptations of two published
adult studies, a realistic adult comparison, a pediatric fixed-helmet versus head-adaptive
extension, and epilepsy (interictal spike) detection and localization examples in both.
The goal and milestones are in `GOAL.md`; the plan, status and decisions log in `PLAN.md`.

This is a proposed study in development. It does not claim that any author of the reproduced
papers has reviewed or approved it. The OPM advantage is tested, not assumed.

## Status

| Milestone | Label | Status |
|---|---|---|
| G0 audit, provenance, plan | - | done |
| G1A Jas et al. 2026 analytical benchmark | REPRO | done, independently reviewed |
| G1B Hunold et al. 2016 depth-orientation spike SNR (MEG part) | ADAPT (+ NEW OPM column) | done, independently reviewed |
| G1C Goldenholz et al. 2009 cortical SNR maps (MEG part) | ADAPT (+ NEW OPM extension) | done, independently reviewed |
| G2 realistic adult OPM vs Neuromag | NEW | done, independently reviewed; frozen as `adult-baseline-v1`, corrected in `adult-baseline-v2` |
| G3 pediatric extension | NEW (size benchmark REPRO) | G3A size benchmark done; G3B done (24-, 18- and 12-month infant templates and two scaled-adult size controls), after the adult freeze `adult-baseline-v2`; the 2-year-template pass independently reviewed (approve with notes; fixes re-verified); the 18- and 12-month templates added on 2026-10-01 |
| G4 epilepsy detection and localization | NEW | adult done; pediatric done (three templates and both size controls); head motion and OPM slippage: bounded secondary extension done; adult and 2-year-template parts independently reviewed |
| G5 software, reproduction, report | - | 124 unit tests pass (also in a clean clone of 57a025e: 123 pass, the full-resolution lead-field check skips without the local cache); `scripts/run_all.sh`; local report in `site/_build` (link-checked, not deployed); private repository, no Pages; clean-environment smoke test passed on 2026-09-30 at a commit before 81168f3 (87 tests then); release prepared, not executed (`docs/release_checklist.md`) |

Adult findings so far (`adult-baseline-v2`), conditional on one adult head (MNE sample subject),
an assumed OPM noise of 15 fT/sqrt(Hz), OPM sensors with no helmet-to-scalp gap beyond the 7-mm
standoff and Neuromag at its
measured (not best) head position; details and caveats in `docs/methods.md` and
`results/g2/G2_report.md`:
- With modelled brain noise (calibrated on gradiometers; it predicts 0.73x the measured
  magnetometer brain noise), a dense on-scalp OPM array (212 single-axis sensors) has 1.13x
  [1.10-1.16] the known-topography detectability of the Neuromag system: about 1.7x for sources
  10-15 mm below the scalp, falling to 1.02-1.05x below 35 mm. An OPM array at Neuromag's own 97
  sites ties it (1.00x), and falls behind after external-field projection (0.90x [0.81-0.95]).
- The advantage disappears with worse assumptions: 1.01x at an OPM noise of 30 fT/sqrt(Hz),
  0.97x with 30 fT/sqrt(Hz) and a 3-mm scalp gap, 0.93x with a 6-mm gap. Without brain noise,
  Neuromag wins (0.76x). The 3- and 1-layer head models agree (1.13x on both).
- Simulated interictal spikes: at 10-20 mm depth the dense array needs about a third less source
  strength for 50 % detection (35 [28-50] vs 53 [44-66] nAm; the intervals overlap, the paired
  strength ratio is 1.51 [1.25-1.66]) and detects more events at 16 of the 18 locations, none
  favouring Neuromag (p < 0.001, uncorrected). Deeper, no location-level difference is established
  (at 20-30 mm, p = 0.23 here but 0.03 in an earlier run whose noise draws differed: not robust).
  The matched array is ahead by location only at 10-20 mm (9/3, p = 0.04 uncorrected; 0.14 in the
  earlier run: not robust) and behind at 45-70 mm (0/7, p = 0.016 uncorrected; the same direction in
  the earlier run). v1's claim that the dense array was favoured in every depth band (14/0, 10/2,
  12/2, 13/0) is superseded. The 72 locations are frontal-heavy (24 frontal, 2 occipital).
- Bounded localization (24 locations, one event each): dipole errors are similar across arrays
  (about 4-5 mm, limited by a 2-mm/2-deg coregistration error); for extended 320-nAm sources the
  distributed (dSPM) estimate tends to be about 5 mm more accurate with either OPM array (p = 0.14
  and 0.046 uncorrected here, 0.003 and 0.01 in an earlier run with different noise draws: the
  direction is stable, the significance is not). Differences not detected are not excluded.
- `adult-baseline-v1` reported 1.21x; that value was inflated by OPM cells reaching into the
  coarse BEM head surface, a numerical error corrected in v2 (see PLAN.md).

Pediatric findings (G3B, NEW), conditional on one adult head, three average infant templates of
one database (24, 18 and 12 months; O'Reilly et al. 2021) and two scaled copies of the adult; the
same Neuromag helmet, sensors and noise for every head; OPM arrays refitted to each head with the
adult rules; details in `docs/methods.md` section 10 and `results/g3b/G3B_report.md`:
- With the child raised to 20-mm contact with the top of the fixed helmet, the dense OPM array's
  known-topography detectability relative to Neuromag combined (D) rises from +0.99 dB in the adult
  to +1.49 dB (school-age-size control), +2.15 dB (2-year-size control), +1.85 dB (2-year
  template), +1.90 dB (18-month template) and +2.06 dB (12-month template): Delta = D_child -
  D_adult = +0.47 [+0.38, +0.58], +1.13 [+0.99, +1.25], +0.73 [+0.44, +1.16], +0.85 [+0.49, +1.27]
  and +1.03 [+0.66, +1.44] dB. Against the gradiometers or magnetometers alone, after external-field
  projection, for the matched-site OPM array and for extended sources the sign is the same.
- The gain comes from the fixed helmet's fit: left at the adult's ear-line position Delta is
  +1.15 to +2.49 dB; laterally centred or at true 18-mm contact it stays +0.46 to +1.16 dB; in a
  counterfactual helmet scaled with the head, centred laterally, it is -0.17, -0.30, -0.37, -0.29
  and -0.12 dB (without the lateral centring the off-centre templates kept +0.28, +0.20 and +0.08
  dB). With the background fixed per unit cortical area, both systems' detectability rises in the
  smaller heads and the on-scalp OPM's rises more (vertex-wise +1.15 and +1.74 dB vs Neuromag +0.65
  and +0.54 dB in the scaled controls); with a helmet that fits, the SQUID gains about as much. The
  children's arrays also have fewer OPM sites (171, 155, 151, 157, 145 vs 212); at an equal site
  count Delta would be larger.
- Delta stays positive for OPM noise 7-30 fT/sqrt(Hz), background variance x0.5 or x2 and a
  1-layer head model; at 30 fT/sqrt(Hz) the adult's D is -0.04 dB (a tie) and the templates' +0.61,
  +0.61 and +0.77 dB. In the scaled controls the gain is largest near the surface (vertex-wise
  +2.25 dB at 10-15 mm for the 2-year size). In the templates the within-stratum gain is largest for
  deep and radial sources, but much of their pooled gain reflects their shallower cortex
  (reweighted to the adult's depth mix, +0.33, +0.54 and +0.32 instead of +0.85, +0.91 and +1.07 dB
  at target level); it is also regionally asymmetric, with the templates off-centre in the helmet.
- At 100 nAm (detectability >= 5, an operational threshold) both systems reach 66 % of the
  adult's and 75-79 % of the templates' usable cortex, the OPM alone a further 2 % and 5 %, and the
  SQUID alone at most 0.15 %.
- Simulated spikes (same detectors and seeds as the adult): the full OPM array's superficial
  detection advantage is about the same in the adult and every smaller head (strength for 50 %
  detection at 10-20 mm, Neuromag combined / dense OPM: adult 53/35; children 52/37, 59/36, 45/30,
  36/27 and 40/29 nAm; paired strength ratio 1.33-1.64 against the adult's 1.51; 14-16 of 18
  locations favour the OPM in each); deeper, a location-level difference appears only for the
  2-year size control (30-70 mm) and the 18-month template (45-70 mm), uncorrected and not
  surviving a correction, so G3B's deeper gain is mostly not resolved at this sample size.
  Localization: dipole errors similar; dSPM of strong focal spikes about 5 mm better with the
  dense OPM in all three templates (p = 0.035, 0.002 and 0.005, uncorrected; 16 comparisons per OPM
  array and anatomy: of these three only the 18-month value survives a Bonferroni correction within
  its anatomy; seven pediatric localization comparisons do in all, and only one, the 12-month dense
  dSPM of strong patches, p = 0.00025, across all children).

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
  about 0.02 deg RMS of in-band rotation in a 1-nT residual field without correction, and at 0.4 deg
  after a homogeneous or 8-term field projection with 1-deg/1-% sensor calibration errors (0.13-0.15
  deg with 3 deg/3 %); the thresholds scale inversely with the field (0.04 deg in 10 nT), and a
  homogeneous projection alone leaves the gradient term (1 dB at 0.14-0.21 deg in a 1-nT/m
  gradient). If the artefact is modelled (oracle covariance), the loss stays below 0.15 dB.
  These are bounds on two mechanisms, not a motion-robustness result; real head-motion statistics,
  sensor dynamic range and gain changes are not modelled.

Labels: REPRO = reproduction with the paper's definitions; ADAPT = adaptation where original
data or details are unavailable; NEW = new experiment or study choice.

## Layout

| Path | Content |
|---|---|
| `src/opmsquid/` | package: sphere model, sensors (Neuromag, OPM), anatomy (adult, infant template, scaled controls), forward models (cached), noise, metrics, background, environment, paper-specific modules (`hunold`, `goldenholz`), G2 comparison, pediatric helmet placements (`pediatric`) |
| `scripts/` | one driver per milestone (`g1a_*`, `g1b_*`, `g1c_*`, `g2_*`, `g3a_*`, `g3b_*`, `g4_*`), full-resolution forward jobs, Fig. 6 digitiser, BEM accuracy studies, `run_all.sh` |
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
bash scripts/run_all.sh     # tests, then every milestone in order (~8 h on a 10-core laptop: adult ~2 h, five pediatric G4 runs ~1 h each)
```

Each driver writes `results/<milestone>/` and records the code commit it ran in its JSON.

The sample data must end up in `data/external/MNE-sample-data` (`src/opmsquid/paths.py`;
override with `OPMSQUID_DATA`). Caches go to `cache/` (`OPMSQUID_CACHE`).

## Local report

`.venv/bin/python scripts/build_site.py` builds a static report into `site/_build/` (git-ignored;
open `site/_build/index.html`) from the committed result files, in the order the goal asks for:
adult benchmarks, realistic adult results, pediatric extension, epilepsy, then methods,
parameters and reproduction. Every number is read from `results/`; the downloads list each file's
size, SHA-256 and the code commit recorded in it; the build fails on any broken link. Nothing is
deployed: repository visibility and GitHub Pages stay unchanged until the owner approves a
release (release-ready and publicly deployed are separate statuses).

## Data and privacy

- The infant templates (24, 18 and 12 months; O'Reilly et al. 2021, from the Neurodevelopmental
  MRI Database of Richards et al. 2016; LGPL-2.1 repository) are not committed. Figures and tables
  derived from them (every `results/g3b` figure that shows a template, `results/g3b/g3b_targets_infant*`,
  `results/g4/*infant*` and the motion results) cite both papers; confirm their redistribution with
  the owner before any public release (the source database has its own terms).
- Not committed: the reference PDFs, the MNE sample data and anatomy (`data/`), computed caches
  (`cache/`). Results contain only derived quantities of the public MNE sample dataset.
- Release preparation (`docs/release_checklist.md`): every figure is set in open fonts (DejaVu
  Sans, STIX); the published Jas et al. Fig. 3 raster used by the legacy verification is CC-BY 4.0
  and carries its attribution; output summaries no longer record local absolute paths. Still the
  owner's decisions: a licence (none yet), the author name and institutional e-mail in the commit
  metadata, the redistribution of template- and fsaverage-derived figures, and two documents that
  name local folders.
- The repository is private. No website is deployed; publication needs explicit owner approval.
