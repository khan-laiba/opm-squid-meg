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
| G3 pediatric extension | NEW (size benchmark REPRO) | not started (needs pediatric anatomy) |
| G4 epilepsy detection and localization | NEW | adult done; pediatric waits for G3 |
| G5 software, reproduction, report | - | in progress (tests, `scripts/run_all.sh`; local report pending) |

Adult findings so far (details and caveats in `docs/methods.md` and `results/g2/G2_report.md`):
with realistic brain noise, a dense on-scalp OPM array (215 single-axis sensors) has about 1.2x
the known-topography detectability of the Neuromag system (about 1.7x for sources 10-15 mm below
the scalp; 1.05x at an OPM noise of 30 fT/sqrt(Hz); 1.13x with a 1-layer head model), while an
OPM array at Neuromag's own 98 sites ties it; without brain noise, Neuromag wins. For simulated
interictal spikes the dense array needs 15-35 % less source strength for 50 % detection,
depending on depth, while the localization accuracy of detected spikes is similar across arrays.

Labels: REPRO = reproduction with the paper's definitions; ADAPT = adaptation where original
data or details are unavailable; NEW = new experiment or study choice.

## Layout

| Path | Content |
|---|---|
| `src/opmsquid/` | package: sphere model, sensors (Neuromag, OPM), anatomy, forward models (cached), noise, metrics, background, environment, paper-specific modules (`hunold`, `goldenholz`), G2 comparison |
| `scripts/` | one driver per milestone (`g1a_*`, `g1b_*`, `g1c_*`, `g2_*`, `g3a_*`, `g4_*`), full-resolution forward jobs, Fig. 6 digitiser, `run_all.sh` |
| `configs/` | paper and study configurations, every value tagged as printed, chosen or digitised |
| `tests/` | `unittest` suite (physics, units, metrics, geometry, noise) |
| `results/` | figures, CSV and JSON per milestone (JSON records the code commit and versions) |
| `docs/` | methods, provenance register, audit, structured literature extractions |
| `legacy/` | earlier Fig. 3 replication work, kept unchanged |

## Reproduce

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -c "import mne; mne.datasets.sample.data_path(path='data/external')"  # ~1.6 GB download
bash scripts/run_all.sh     # tests, then every milestone in order (~2-3 h; lead fields ~6 min each)
```

Each driver writes `results/<milestone>/` and records the code commit it ran in its JSON.

The sample data must end up in `data/external/MNE-sample-data` (`src/opmsquid/paths.py`;
override with `OPMSQUID_DATA`). Caches go to `cache/` (`OPMSQUID_CACHE`).

## Data and privacy

- Not committed: the reference PDFs, the MNE sample data and anatomy (`data/`), computed caches
  (`cache/`). Results contain only derived quantities of the public MNE sample dataset.
- Some G1A figure PDFs embed licensed font subsets to match the original artwork
  (`results/g1a/fig3/README.md`); regenerate with open fonts before any public release.
- The repository is private. No website is deployed; publication needs explicit owner approval.
