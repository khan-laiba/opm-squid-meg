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
| G3 pediatric extension | NEW (size benchmark REPRO) | G3A size benchmark done; G3B done (2-year infant template and two scaled-adult size controls), after the adult freeze `adult-baseline-v2` |
| G4 epilepsy detection and localization | NEW | adult done; pediatric done (2-year template, school-age size control) |
| G5 software, reproduction, report | - | in progress (tests, `scripts/run_all.sh`, local report; clean-environment smoke test passed on 2026-09-30 at a commit before 81168f3, with the 87 tests of the time) |

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

Pediatric findings (G3B, NEW), conditional on one adult head, one 2-year average template
(O'Reilly et al. 2021) and two scaled copies of the adult; the same Neuromag helmet, sensors and
noise for every head; OPM arrays refitted to each head with the adult rules; details in
`docs/methods.md` section 10 and `results/g3b/G3B_report.md`:
- With the child raised to 20-mm contact with the top of the fixed helmet, the dense OPM array's
  known-topography detectability relative to Neuromag combined (D) rises from +0.99 dB in the adult
  to +1.49 dB (school-age-size control), +2.15 dB (2-year-size control) and +1.85 dB (2-year
  template): Delta = D_child - D_adult = +0.47 [+0.38, +0.58], +1.13 [+0.99, +1.25] and +0.73
  [+0.44, +1.16] dB. Against the gradiometers or magnetometers alone, after external-field
  projection, for the matched-site OPM array and for extended sources the sign is the same.
- The gain comes mainly from the helmet's fit: left at the adult's ear-line position Delta is
  +1.15 to +1.79 dB; in a counterfactual helmet scaled with the head it is -0.17, -0.30 and
  +0.28 dB. In the fixed helmet the child's cortex is farther from the SQUIDs, whose brain noise
  falls towards their intrinsic floor while the on-scalp OPM's does not.
- Delta stays positive for OPM noise 7-30 fT/sqrt(Hz), background variance x0.5 or x2 and a
  1-layer head model; at 30 fT/sqrt(Hz) the adult's D is -0.04 dB (a tie) and the template's
  +0.61 dB. In the template the gain is concentrated deeper than 30 mm and is regionally asymmetric
  (it follows the head's offset in the helmet).
- At 100 nAm (detectability >= 5, an operational threshold) both systems reach 66 % of the
  adult's and 75 % of the template's usable cortex, the OPM alone a further 2 % and 5 %, and the
  SQUID alone none.
- Simulated spikes (same detectors and seeds as the adult): the full OPM array's superficial
  detection advantage is about the same in the adult and both smaller heads (strength for 50 %
  detection at 10-20 mm, Neuromag combined / dense OPM: 53/35, 52/37 and 45/30 nAm; 16 of 18
  locations favour the OPM in each, and 0, 0 and 2 favour Neuromag); deeper, no location-level
  difference is established for the
  practical detector, so G3B's deeper gain is not resolved at this sample size. Localization: dipole
  errors similar; dSPM of strong focal spikes in the template 7.5 vs 14.8 mm (p = 0.035,
  uncorrected).

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
| `legacy/` | earlier Fig. 3 replication work, kept unchanged |

## Reproduce

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -c "import mne; mne.datasets.sample.data_path(path='data/external')"  # ~1.6 GB download
.venv/bin/python -c "import mne; mne.datasets.fetch_infant_template('2yr', subjects_dir='data/external/infant_subjects')"  # ~392 MB (G3B, pediatric G4)
bash scripts/run_all.sh     # tests, then every milestone in order (~2 h, of which ~40 min for lead fields)
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

- The 2-year infant template (O'Reilly et al. 2021, from the Neurodevelopmental MRI Database of
  Richards et al. 2016; LGPL-2.1 repository) is not committed. Figures derived from it
  (`results/g3b/*template*`, `*infant2yr*`) cite both papers; confirm their redistribution with the
  owner before any public release (the source database has its own terms).
- Not committed: the reference PDFs, the MNE sample data and anatomy (`data/`), computed caches
  (`cache/`). Results contain only derived quantities of the public MNE sample dataset.
- Before any public release: PDFs in `results/g1a/` and `legacy/` embed licensed font subsets
  (Myriad Pro, Times New Roman) to match the original artwork (`results/g1a/fig3/README.md`), so
  regenerate them with open fonts; the published Jas et al. Fig. 3 raster used by the legacy
  replica is CC-BY 4.0 and needs attribution; commit metadata carries the author's name and
  institutional e-mail, and a few tracked files contain local absolute paths.
- The repository is private. No website is deployed; publication needs explicit owner approval.
