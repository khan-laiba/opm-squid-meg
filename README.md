# OPM vs SQUID MEG: adult-to-pediatric simulation study

Simulation study comparing on-scalp optically pumped magnetometers (OPM) with the Neuromag
(Vectorview/TRIUX-type) SQUID system: an analytical benchmark, adaptations of two published
adult studies, a realistic adult comparison, and (planned) pediatric and epilepsy extensions.
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
| G2 realistic adult OPM vs Neuromag | NEW | done, independently reviewed; frozen as `adult-baseline-v1` |
| G3 pediatric extension | NEW (size benchmark REPRO) | G3A size benchmark done; G3B needs pediatric anatomy |
| G4 epilepsy detection and localization | NEW | adult done; pediatric waits for G3 |
| G5 software, reproduction, report | - | in progress (tests, `scripts/run_all.sh`, local report; clean-environment smoke test pending) |

Adult findings so far, conditional on one adult head (MNE sample subject), an assumed OPM noise
of 15 fT/sqrt(Hz), OPM sensors with no gap to the scalp and Neuromag at its measured (not best)
head position; details and caveats in `docs/methods.md` and `results/g2/G2_report.md`:
- With modelled brain noise (calibrated on gradiometers; it predicts 0.73x the measured
  magnetometer brain noise), a dense on-scalp OPM array (215 single-axis sensors) has 1.21x
  [1.19-1.23] the known-topography detectability of the Neuromag system, about 1.7x for sources
  10-15 mm below the scalp. An OPM array at Neuromag's own 98 sites ties it (1.02x), and falls
  behind after external-field projection (0.92x [0.85-0.95]).
- The advantage shrinks or disappears with worse assumptions: 1.05x at an OPM noise of
  30 fT/sqrt(Hz), 1.00x with 30 fT/sqrt(Hz) and a 3-mm scalp gap, 0.95x with a 6-mm gap, and
  1.13x with a 1-layer head model. Without brain noise, Neuromag wins (0.77x).
- Simulated interictal spikes: the dense array's point estimates of the strength for 50 %
  detection are 15-35 % lower, depending on depth, with wide intervals (for example 57 [39-113]
  vs 67 [50-145] nAm at 20-30 mm). By location, the practical detector favours the dense array
  in every depth band; the matched array does not differ.
- Bounded localization (24 locations, one event each): no consistent difference between arrays
  was detected for detected spikes; differences are not excluded.

Labels: REPRO = reproduction with the paper's definitions; ADAPT = adaptation where original
data or details are unavailable; NEW = new experiment or study choice.

## Layout

| Path | Content |
|---|---|
| `src/opmsquid/` | package: sphere model, sensors (Neuromag, OPM), anatomy, forward models (cached), noise, metrics, background, environment, paper-specific modules (`hunold`, `goldenholz`), G2 comparison |
| `scripts/` | one driver per milestone (`g1a_*`, `g1b_*`, `g1c_*`, `g2_*`, `g3a_*`, `g4_*`), full-resolution forward jobs, Fig. 6 digitiser, `run_all.sh` |
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

- Not committed: the reference PDFs, the MNE sample data and anatomy (`data/`), computed caches
  (`cache/`). Results contain only derived quantities of the public MNE sample dataset.
- Before any public release: PDFs in `results/g1a/` and `legacy/` embed licensed font subsets
  (Myriad Pro, Times New Roman) to match the original artwork (`results/g1a/fig3/README.md`), so
  regenerate them with open fonts; the published Jas et al. Fig. 3 raster used by the legacy
  replica is CC-BY 4.0 and needs attribution; commit metadata carries the author's name and
  institutional e-mail, and a few tracked files contain local absolute paths.
- The repository is private. No website is deployed; publication needs explicit owner approval.
