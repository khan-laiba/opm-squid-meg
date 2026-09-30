# Parameter, provenance and assumptions register

Every numerical choice used by the study, with its origin. Categories never mix:

- **J / HU / GO**: reported in Jas et al. 2026 / Hunold et al. 2016 / Goldenholz et al. 2009
  (page numbers refer to `docs/literature/*.md`).
- **HW**: hardware specification or vendor/site documentation.
- **R**: recovered by measurement from a published figure, not stated in the text.
- **A**: new study assumption (with rationale and the sensitivity analysis that covers it).
- **D**: study decision (design choice, not a physical parameter).
- **U**: unresolved ambiguity (kept open, not silently repaired).

Status: `verified` (checked against the source or reproduced numerically), `transcribed`
(copied from GOAL.md; primary source not available here), `assumed`, `open`.

## Jas et al. 2026 (analytical benchmark)

| ID | Parameter | Value | Source | Status | Used in |
|---|---|---|---|---|---|
| J-h | Adult head radius h | 95 mm | GOAL.md; paper Table 1 "Adult" | verified (legacy replica) | G1A |
| J-b | Adult brain radius b | 80 mm | GOAL.md; paper Table 1 | verified | G1A |
| J-Q | Tangential dipole moment | 30 nAm | paper | verified | G1A |
| J-xi | Scalp standoff OPM / SQUID | 0 / 18 mm | paper | verified | G1A |
| J-eta | Relative RMS noise ratio sigma_OPM / sigma_SQUID | 1-6 (3 in Fig. 3) | paper | verified | G1A |
| J-eq1 | Peak radial field of a tangential dipole (Eq. 1) | closed form, see `opmsquid.sphere.bmax_radial` | paper Eq. 1 | verified vs independent Sarvas maximum (< 1e-8) | G1A |
| J-deq | Equal-SNR depth at eta = 3 | 27.665 mm (paper: "approximately 28 mm") | Eq. 3, computed | verified | G1A |
| J-opm-standoff | OPM sensing-centre distance from the helmet inner surface | ~7 mm (2-mm helmet shell + half of a 10-mm cell) | paper p. 10 | verified (extraction) | G2 (A-OPM-STANDOFF) |
| J-opm-cell | OPM vapour-cell size | 10-mm cube, FieldLine Gen2, single axis, normal component | paper p. 10 | verified (extraction) | G2 (A-OPM-CELL) |
| J-table1 | Head models (h, b) | newborn (55, 48), 1 y (70, 62), 8 y (85, 73), adult (95, 80) mm | paper Table 1, p. 7 | verified (extraction) | G3 size benchmark |
| J-noise-lit | Literature intrinsic noise | SQUID 2-5 fT/sqrt(Hz), OPM 7-30 fT/sqrt(Hz) (Brookes et al. 2022) | paper p. 4 | verified (extraction) | A-OPM-NOISE sweep range |
| J-eta-meas | Measured noise ratio (not stated) | ~4.6 (Fig. 8B, sigma_OPM / sigma_SQUID at the closest runs) | derived from the figure | derived | context |
| J-toy | Toy experiment (Fig. 6) | target r_Q = 0.6 b; noise dipoles 0.4 b and 0.8 b; all 30 nAm +y; SNR = ratio of peak fields; xi 0-60 mm | pp. 9, 16-17 | reproduced (G1A) | G1A |

## Recovered from the published Fig. 3 raster (legacy replica)

| ID | Parameter | Value | Status |
|---|---|---|---|
| R-sigma | Absolute SQUID noise sigma_SQUID | B_max(0.8 b, h + 18 mm) = 0.35458 pT (SNR_SQUID = 1 at d = 31 mm); sigma_OPM = eta sigma_SQUID | recovered; the paper gives only eta. Absolute SNR values depend on it, ratios and d_eq do not. |
| R-grid | Curve sampling | r_Q = linspace(0, b, 101)[1:] (d = 15 ... 94.2 mm) | recovered |
| R-marker | Dotted d_eq line position | 27.530 mm (deepest point with SNR_OPM > SNR_SQUID on d = linspace(15, 95, 250)) vs exact root 27.665 mm | recovered; authors' exact rule not identifiable (U-J1) |

## Hardware

| ID | Parameter | Value | Source | Status | Used in |
|---|---|---|---|---|---|
| HW-mag-noise | Neuromag TRIUX magnetometer typical white noise | 3.5 fT/sqrt(Hz) = 3.5e-15 T/sqrt(Hz) | TRIUX specification image (user-supplied) | transcribed from GOAL.md; image not on disk (U-HW1) | G2 |
| HW-grad-noise | Neuromag TRIUX planar gradiometer typical white noise | 3.6 fT/(cm sqrt(Hz)) = 3.6e-13 T/(m sqrt(Hz)) | same | transcribed (U-HW1) | G2 |
| HW-18mm | Pickup coil to room-temperature (helmet) surface | 18 mm | same | transcribed. Not a scalp gap: measured scalp-to-coil distances for the sample subject are 25-41 mm (median 31 mm) | G2 |
| HW-T3 | MRN coil types | planar gradiometer 3014, magnetometer 3024 (T3) | mrn.org (read 2026-09-30) | verified | G2 |
| HW-MRN-floor | MRN shielded-room noise floor | "5-7 fT" (bandwidth unstated) | mrn.org | transcribed; not used as a spectral density | context |
| HW-geometry | Neuromag sensor geometry | MNE sample recording (MGH Vectorview, 306 channels; stored coils 3012 x 204, 3024 x 102) | MNE 1.13.2 sample dataset | verified (read) | G1B, G1C, G2 |
| HW-coildef | Coil integration | MNE 1.13.2 `coil_def.dat`, 'accurate' level (3014: 8 points, 26.39-mm coils, 16.80-mm baseline; 3024: 16 points, 21-mm coil) | MNE | verified | all |

## New study assumptions (A) and decisions (D)

| ID | Item | Value / rule | Rationale and coverage |
|---|---|---|---|
| D-T3 | Gradiometer coil variant | 3012 in the sample file replaced by 3014 | Same geometry and integration points in MNE 1.13.2 except the 1-point model (0.3 mm offset); matches MRN (HW-T3). |
| A-OPM-CELL | OPM sensing volume | 10-mm cube, custom coil 9901 (27-point Gauss at MNE 'accurate') | J-opm-cell; MNE forward = independent quadrature (< 1e-6). Effect at the field peak: -1.9 % (10 mm source distance), -0.4 % (15 mm), -0.04 % (30 mm). |
| A-OPM-STANDOFF | Sensing centre from helmet inner surface | 7 mm | J-opm-standoff |
| A-OPM-GAP | Helmet inner surface to MRI scalp | 0 mm baseline; swept (e.g. 0, 3, 6 mm) | Hair and helmet fit are not modelled otherwise. |
| A-OPM-AXIS | OPM sensitive axis | smoothed scalp normal (10-mm radius), single axis | Multi-axis OPMs are a secondary extension. |
| A-OPM-COVER | OPM coverage | scalp above the plane through LPA, RPA and nasion + 3 cm | Excludes face and neck; sample subject: 99 of 102 matched sites kept. |
| A-OPM-CLEAR | Physical clearance | sensing centre >= standoff - 1 mm from every scalp point, sites moved outward if needed (sample: 3 sites, +3 to +4 mm over the ear pinnae and brow) | Package collisions with the pinna/brow. |
| A-OPM-NOISE | OPM intrinsic noise | sweep 7-30 fT/sqrt(Hz) (white) | No single verified device specification (GOAL.md). |
| D-ANAT | Primary adult anatomy | MNE `sample` subject | Individual adult MRI with real helmet position; fsaverage secondary. |

## Unresolved ambiguities (U)

| ID | Item | Handling |
|---|---|---|
| U-J1 | Jas Fig. 3 d_eq marker at 27.53 mm vs Eq. 3 root 27.665 mm vs "28 mm" in the text | Both values reported; exact root used in analyses. |
| U-J2 | Fig. 4B d_eq printed 34 mm, drawn ~35.0 mm, exact 35.67 mm; Fig. 4C printed 19 mm, drawn ~19.8, exact 19.82 | Exact roots used; printed and drawn values listed in `results/g1a/g1a_benchmark.json`. |
| U-J3 | eta0 / eta1 printed 1.7 / 5.3; exact 1.6829 / 5.3086 (a deep crossing exists for 1.683 < eta < 1.7) | Exact values used. |
| U-J4 | Fig. 5B "normalized d_eq" = (d_eq - (h - b))/b (depth below the brain surface as a fraction of b), not d_eq/b; text speaks of brain *volume* | G3 reproduces the radial fraction as plotted and reports the true volume fraction separately (adult 40.4 %, newborn 87.2 % at eta = 3). |
| U-J5 | Fig. 5 uses s = h + 18 mm for every head (size-following shell) | Reproduced as the "size-following, constant-standoff benchmark"; the fixed-helmet experiment is separate (G3). |
| U-J6 | Fig. 6 depths: text 32/48/64 mm (these are r_Q = 0.4b/0.6b/0.8b) vs caption 63/47/31 mm; noise-dipole moments unstated (30 nAm reproduces the figure) | Caption values and 30 nAm used. |
| U-J7 | sigma_N20 defined as trial SD (p. 12) vs "standard error of mean" (Fig. 8) | Experimental context only (no recordings available). |
| U-J8 | Experimental OPM ~7 mm from the scalp (2-mm shell + half of a 10-mm FieldLine Gen2 cell) vs simulated OPM at xi = 0 | G1A keeps xi = 0 (paper model); G2 uses 7 mm (A-OPM-STANDOFF). |
| U-HW1 | TRIUX specification image not available on disk | Values transcribed in GOAL.md used and flagged; to verify when the image is supplied. |
