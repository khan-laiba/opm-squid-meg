# Sensor noise, array design and helmet fit in on-scalp and cryogenic MEG

This repository contains the code, the stored analysis outputs and the manuscript sources for

> **Sensor noise, array design and helmet fit in on-scalp and cryogenic magnetoencephalography: a simulation study
> in adult and pediatric head models.**
> Laiba Khan and Mainak Jas. Manuscript prepared for submission to *NeuroImage*, 2026.

**Paper:** [HTML](https://khan-laiba.github.io/opm-squid-meg/) · [PDF](paper/manuscript.pdf) ·
[Supplementary material](paper/supplementary.pdf)

## Overview

Optically pumped magnetometers (OPMs) can be placed on the scalp, where the magnetic field of cortical sources is
stronger than at the superconducting quantum interference devices (SQUIDs) of conventional magnetoencephalography
(MEG). Whether the stronger field yields a higher signal-to-noise ratio (SNR) depends on the noise: the white noise of
OPMs is higher than that of SQUID magnetometers (15 against 3.5 fT/√Hz in the primary model), and the field of
spontaneous brain activity, the brain noise, is also stronger at the scalp. A spherical-head analysis with point
sensors at fixed standoffs and an assumed ratio of OPM to SQUID noise
([Jas et al., 2026](https://doi.org/10.64898/2026.08.17.744953)) predicted the largest on-scalp benefit for superficial
sources and small heads.

This study examines that prediction in realistic head models. The 306-channel Neuromag SQUID system is compared with
two arrays of single-axis OPMs, a dense array (208 sensors on the adult head) and an array at 98 of the 102 Neuromag
sensor sites (Fig. 1), by the detectability of each cortical source: the output SNR of a matched filter that knows the
field pattern of the source and the covariance of the sensor, brain and environmental noise. Forward fields are
computed with three-layer boundary-element models of an adult head (the MNE sample subject), two copies of it scaled
down to school-age size and to the head circumference of the 24-month template, and the 24- and 12-month infant
templates;
the 18-month template and three school-aged children are reported separately (Section S10). Each head is placed in the
adult SQUID helmet and in a copy of that helmet scaled until its median magnetometer-to-scalp distance equals the
adult's, 29.6 mm (Fig. 2). The modeled Neuromag noise is compared with the measured noise covariance of the sample
recording, and interictal spikes are simulated and detected by idealized matched-filter detectors.

<p align="center">
  <img src="https://khan-laiba.github.io/opm-squid-meg/images/Figure_1.png"
       alt="Fig. 1: the Neuromag sensor sites and the two OPM arrays on the adult head" width="49%">
  <img src="https://khan-laiba.github.io/opm-squid-meg/images/Figure_2.png"
       alt="Fig. 2: the fixed adult helmet and the helmet fitted at the adult's gap" width="49%">
</p>

*Left, Fig. 1: the 102 Neuromag sensor sites and the two OPM arrays on the adult head. Right, Fig. 2: sagittal and
coronal sections through the adult, the 12-month infant template and child B (8.3 years; Section S10), with the
Neuromag magnetometers of the fixed adult helmet at top contact (black), those of the helmet fitted at the adult's
gap (vermilion) and the dense OPM sites (blue).*

## Main results

In items 1–3, ratios are medians over the adult's 7,661 cortical sources of the source-wise detectability ratio
OPM/Neuromag, with 95 % parcel-bootstrap intervals. In item 4, *D* is that ratio in decibels and Δ the paired change
in *D* from the adult (area-weighted medians).

1. With sensor noise alone, Neuromag has the higher detectability in the adult (ratio 0.74 for the dense array, 0.53
   for the site-matched array). Brain noise, the dominant noise in both systems, reverses this ranking for the dense
   array: with sensor and brain noise and an OPM white noise of 15 fT/√Hz, its ratio is 1.14 (1.12–1.17), whereas
   that of the site-matched array is 1.01 (0.99–1.03). Fig. 3, Table 3.
2. The OPM noise level and the array design determine which system has the higher detectability: under the modeled
   Neuromag noise, the ratio reaches 1 at an OPM noise of 31.5 (28.5–34.9) fT/√Hz for the dense array and 16.7
   (14.0–19.4) fT/√Hz for the site-matched array. The Neuromag channels used for comparison and the SNR measure
   determine which system a sensor-level SNR favors. Figs. 3 and 4D, Table 3.
3. The dense-array advantage is as large as the error of the noise model for Neuromag: with its measured noise
   covariance, Neuromag has 1.14 (1.08–1.20) times its modeled detectability. Because the OPM covariance was not
   measured, the implied ratio was computed under five assumptions about how each system's noise departs from the
   model. It ranges from 0.80 to 1.14 for the dense array and is 1.03 (0.99–1.08) with Neuromag's noise as measured
   and the OPM noise as modeled, in which case the dense array's break-even level would fall to near 15 fT/√Hz.
   Fig. 4, Table 3.
4. In the fixed adult helmet, the four principal smaller heads have larger median magnetometer-to-scalp distances
   than the adult (34.1–40.4 against 28.4 mm) and a larger dense-array advantage (Δ = +0.44 to +1.07 dB). Relative to
   the fitted helmet, the fixed helmet lowers Neuromag's detectability more in each of these heads than in the adult
   (the interaction, +0.40 to +1.10 dB; independent of the assumed OPM noise). In the fitted helmet, with a dense
   array built on each head, Δ is +0.02 to +0.05 dB; with the adult's array reduced to each head's number of sites,
   it is +0.20 to +0.49 dB. Figs. 8 and 9, Table 4.
5. In the pre-specified spike run, under the modeled noise and at a nominal rate of one false event per minute, the
   dipole strength detected in 50 % of events (S50) for focal spikes 10–20 mm below the scalp is 22–32 % lower with
   the dense array than with Neuromag, in the adult and in each principal smaller head (S50 ratios 1.28–1.47, every
   interval above 1). The detectors are idealized, and the smaller heads were simulated only at top contact in the
   fixed adult helmet. Fig. 10, Table 5.

These are results of simulations. The OPM noise covariance was not measured, and the noise model was calibrated and
checked on a single recording; Section 4.4 of the paper discusses these and the other limitations.

## Repository structure

```
src/opmsquid/   Python package: sensor arrays and helmets, head and forward models, noise model,
                detectability, spike simulation, detection and localization
scripts/        analysis, figure and fact scripts, and run_all.sh (scripts/README.md)
configs/        analysis parameters, the declaration of the pre-specified spike run and the SHA-256
                manifests of the downloaded MRI files
results/        stored outputs of every analysis: JSON and CSV files, each stamped with the commit that
                produced it, and the analyses' own figures and reports
paper/          manuscript and supplement sources, print figures, build scripts and compiled PDFs
docs/           detailed methods, the provenance of every parameter, records of the literature used, and
                commit_map.tsv (commit hashes from before the history rewrite of 4 October 2026 and their
                published equivalents)
tests/          unit tests of the numerical components, the fact scripts and the data fetcher
legacy/         replica of Fig. 3 of Jas et al. (2026), imported by the spherical benchmark of Section S2
```

As declared in the paper (Section 2.7), the analysis code was written and run with the assistance of a generative AI
tool (Claude, Anthropic); its numerical components are covered by the unit tests.

## System requirements

The analyses were run on macOS (Apple M4, 16 GB of memory) with Python 3.12.14, MNE-Python 1.13.2, NumPy 2.5.3 and
SciPy 1.18.1; `requirements.txt` pins the full environment, and other platforms have not been tested. No GPU is
needed. The external data occupy about 5 GB, and the caches written to `cache/` reached 29 GB; the environment
variables `OPMSQUID_DATA` and `OPMSQUID_CACHE` move them. Building the PDFs requires
[tectonic](https://tectonic-typesetting.github.io/) and `pdftotext` (poppler); the HTML version also requires
`pdftoppm`.

## Installation

```bash
git clone https://github.com/khan-laiba/opm-squid-meg.git
cd opm-squid-meg
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Installation takes a few minutes.

## External data

No data are redistributed in this repository. The analyses use three public datasets, downloaded into
`data/external/`:

| Dataset | Source | Terms |
|---|---|---|
| MNE sample dataset (Gramfort et al., 2013) | the release fetched by MNE-Python 1.13.2 (OSF file 86qa2, version 6) | provided for learning the MNE software, not for evaluating the MEG or MRI system |
| Infant templates at 24, 18 and 12 months (O'Reilly et al., 2021; Richards et al., 2016) | MNE-Python's template fetcher | LGPL-2.1 |
| Three school-aged children of OpenNeuro ds005234 v2.2.0 (Fadeev et al., 2024, 2025) | [doi:10.18112/openneuro.ds005234.v2.2.0](https://doi.org/10.18112/openneuro.ds005234.v2.2.0) | CC0 in the dataset's metadata; its authors ask that it be acknowledged as distributed under CC BY with its DOI |

```bash
# MNE sample dataset
.venv/bin/python -c "import mne; mne.datasets.sample.data_path(path='data/external')"
# infant templates
mkdir -p data/external/infant_subjects
.venv/bin/python -c "import mne; [mne.datasets.fetch_infant_template(a, subjects_dir='data/external/infant_subjects') for a in ('2yr', '18mo', '12mo')]"
# school-aged children, checked against SHA-256 sums
.venv/bin/python scripts/fetch_school_subjects.py
```

The MEG recording of the sample dataset is used only to calibrate and check the noise model; no value is offered as
an evaluation of the Neuromag system or of the MRI scanner.

## Reproducing the paper

**From the stored results.** Every result in the manuscript and the supplement is read from `results/` by
`scripts/report_facts*.py`, with parameters from `configs/` and `docs/`; none is typed by hand. The print figures are
committed in `paper/figures/`, so the PDFs can be rebuilt without the external data:

```bash
.venv/bin/python paper/build_paper.py      # paper/manuscript.pdf and paper/supplementary.pdf
```

The build also writes `paper/build/numbers_used.tsv` (not committed), which lists every value with the file and key
it was read from. `.venv/bin/python paper/export_figures.py` redraws the print figures from `results/`; it also reads
the external data (the cortical surfaces for the maps and the T1-weighted images for Fig. S6) and stops without them.
The results reported in the paper are those of commit f167fd3, and `results/` has not changed since.
[`scripts/README.md`](scripts/README.md) gives the scripts and result files behind each figure and table.

**From the external data.** `scripts/run_all.sh` runs the unit tests, then the analyses in the order in which they
depend on one another, and finally the figure and paper builds; `RESUME=1 scripts/run_all.sh` continues an
interrupted run on an unmodified checkout.

```bash
scripts/run_all.sh
```

The stored results were produced by running these steps individually; `run_all.sh` chains them in the same order but
has not been run end to end in one pass. On the machine above, most steps took minutes to an hour. The longest were
the eight pediatric spike runs (40–75 min each) and the pre-specified spike run (about 3 h, 17–22 min per anatomy,
with up to 8.6 GB of memory); the whole pipeline takes the better part of a day. The random seeds are fixed and listed
in Table S17. Every JSON and CSV result file records the commit of the code that produced it; files written before the
history was rewritten on 4 October 2026 (author identity only; contents and dates unchanged) record the earlier
hashes, which [`docs/commit_map.tsv`](docs/commit_map.tsv) maps to the published commits.

**Pre-specified analysis.** The spike run of Sections 2.5 and S7 was declared in
[`configs/g4_confirmatory.toml`](configs/g4_confirmatory.toml), with its endpoint, its test (a two-sided sign-flip test
on the per-location differences in detection counts, Holm-adjusted over the nine anatomies at α = 0.05), its seeds and
its sample size (36 locations per anatomy, five noise replicates), in commit 10e37b9 of 4 October 2026. It was run at
commit ae458a9, which differs from 10e37b9 only in the result files of other analyses. The number of locations was
chosen after two pilot runs in child A whose outputs were not kept. Only the history of this repository attests the
declaration; the run was not registered externally.

**Tests.** After the external data have been downloaded, `.venv/bin/python -m unittest discover -s tests -t .` runs
299 tests (about 10 min). Without the MNE sample dataset, the run reports six errors and skips 35 tests.

## Citation

If you use this code or these results, please cite the paper ([`CITATION.cff`](CITATION.cff)):

```bibtex
@unpublished{khan2026opmsquid,
  author = {Khan, Laiba and Jas, Mainak},
  title  = {Sensor noise, array design and helmet fit in on-scalp and cryogenic magnetoencephalography:
            a simulation study in adult and pediatric head models},
  note   = {Manuscript prepared for submission to NeuroImage},
  year   = {2026}
}
```

## License

The code is released under the MIT license ([`LICENSE`](LICENSE)). The manuscript text and figures, the
documentation and the result files are released under CC BY 4.0, and third-party material keeps its own terms
([`LICENSE-CONTENT.md`](LICENSE-CONTENT.md)).

## Contact

Correspondence about the paper: Mainak Jas (mjas@mgh.harvard.edu). Questions about the code can be raised as
[issues](https://github.com/khan-laiba/opm-squid-meg/issues).
