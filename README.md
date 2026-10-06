# Sensor noise, array design and helmet fit in on-scalp and cryogenic MEG

Code, stored results and manuscript sources for

> Khan L, Jas M. Sensor noise, array design and helmet fit in on-scalp and cryogenic magnetoencephalography: a
> simulation study in adult and pediatric head models. Manuscript prepared for submission to *NeuroImage*, 2026.

**Read the paper:** [HTML](https://khan-laiba.github.io/opm-squid-meg/) ·
[PDF](paper/manuscript.pdf) · [Supplementary material](paper/supplementary.pdf)

## About

The study compares on-scalp optically pumped magnetometers (OPMs) with the 306-channel Neuromag SQUID system in
simulation. A dense OPM array and an array of OPMs at the Neuromag sensor sites are compared with Neuromag through the
matched-filter detectability of cortical sources, with explicit sensor, brain and room noise, in an adult head and in
smaller heads placed in the adult helmet and in a helmet fitted at the adult's gap. The noise model is checked against
the measured noise of the Neuromag recording of the MNE sample dataset, and simulated interictal spikes are detected by
idealized matched-filter detectors.

## Repository contents

| Folder | Contents |
|---|---|
| `src/opmsquid/` | Python package: sensor arrays, head models, forward and noise models, detectability, spike simulation, detection and localization |
| `scripts/` | analysis and figure scripts, and `run_all.sh`, which runs them in order; see [`scripts/README.md`](scripts/README.md) |
| `configs/` | analysis parameters |
| `results/` | stored outputs of every analysis; every number in the paper is read from these files |
| `paper/` | manuscript and supplement sources, print figures, build scripts and compiled PDFs |
| `docs/` | detailed methods notes, the provenance of every parameter, and records of the literature used |
| `tests/` | unit tests |
| `legacy/` | replica of Fig. 3 of Jas et al. (2026), used by the spherical benchmark |

## Installation

The analyses were run with Python 3.12 on macOS.

```bash
git clone https://github.com/khan-laiba/opm-squid-meg.git
cd opm-squid-meg
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Building the PDFs also needs [tectonic](https://tectonic-typesetting.github.io/) and poppler (`pdftotext`, `pdftoppm`).

## Data

No data are redistributed. The analyses use three public datasets, downloaded into `data/external/`:

```bash
# MNE sample dataset (Gramfort et al., 2013)
.venv/bin/python -c "import mne; mne.datasets.sample.data_path(path='data/external')"
# infant templates at 24, 18 and 12 months (O'Reilly et al., 2021)
mkdir -p data/external/infant_subjects
.venv/bin/python -c "import mne; [mne.datasets.fetch_infant_template(a, subjects_dir='data/external/infant_subjects') for a in ('2yr', '18mo', '12mo')]"
# three school-aged children of OpenNeuro ds005234 v2.2.0 (Fadeev et al., 2024, 2025), checked against SHA-256 sums
.venv/bin/python scripts/fetch_school_subjects.py
```

## Reproducing the paper

The figures, tables and every number of the paper are produced from the stored results in `results/`, without
rerunning the simulations:

```bash
.venv/bin/python paper/export_figures.py   # print figures in paper/figures/ (the cortical maps need the data)
.venv/bin/python paper/build_paper.py      # paper/manuscript.pdf and paper/supplementary.pdf
```

`paper/build/numbers_used.tsv` then lists each number in the paper with the result file and key it was read from.
[`scripts/README.md`](scripts/README.md) gives the script and result files behind each figure and table.

To rerun all analyses from scratch (roughly 11 h on a laptop for the main analyses, plus the pre-specified spike runs):

```bash
scripts/run_all.sh
```

The unit tests run with `.venv/bin/python -m unittest discover -s tests -t .`.

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

The code is released under the MIT license ([`LICENSE`](LICENSE)). The manuscript text and figures, the documentation
and the result files are released under CC BY 4.0; third-party material keeps its own terms
([`LICENSE-CONTENT.md`](LICENSE-CONTENT.md)).
