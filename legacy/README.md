# Replica of Fig. 3 of Jas et al. (2026)

`replicate_figure3.py` reproduces Fig. 3 of Jas et al. (2026, bioRxiv, doi:10.64898/2026.08.17.744953): radial point
magnetometers 0 and 18 mm above a 95/80-mm spherical head, a 30-nAm tangential dipole and an OPM noise three times the
SQUID noise. The spherical benchmark of the paper (Section S2, Fig. S1) imports it as a module
(`scripts/g1a_jas_benchmark.py`), and the parameters it recovers from the published figure are documented in its
docstring and in `docs/provenance_register.md`.

`verify_figure3.py` compares the replica's output pixel by pixel with the figure in the preprint PDF, which is not
redistributed here (pass its path with `--pdf`).

Run from this folder with the project environment, e.g. `../.venv/bin/python replicate_figure3.py`; the outputs are
written to `figure3_output/` next to the script.
