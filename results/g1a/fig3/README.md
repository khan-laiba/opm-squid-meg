# Fig. 3 replica (G1A, REPRO)

Reproduction of Jas et al. 2026 Fig. 3 (analytical sphere, eta = 3) by `legacy/replicate_figure3.py`, called from
`scripts/g1a_jas_benchmark.py`. Status and parameters: `../g1a_benchmark.json`.

- `figure3_data.csv`: curve samples (depth, peak fields, SNRs); `figure3_summary.json`: checks against the published raster.
- The dotted d_eq line is drawn at 27.53 mm by a grid rule fitted to the published raster (U-J1); the exact root is 27.665 mm.
- The PNG carries a sub-pixel registration (up to 1.4 px) applied only to overlay the published raster.
- The PDF embeds subsets of Myriad Pro and Times New Roman (licensed fonts, used to match the paper's artwork); the
  SVG carries outlines. Regenerate with open fonts before any public release.
