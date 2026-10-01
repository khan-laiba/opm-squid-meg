# G3B: fixed adult Neuromag helmet versus head-adaptive OPM on smaller heads (NEW)

Code commit: see `g3b_summary.json` (provenance). Configuration: `configs/g3b_pediatric.toml`. Primary placement: `top` (head raised to 20-mm contact). Metric: known-topography detectability in dB; D = OPM - SQUID; Delta = D_child - D_adult. Intervals: parcel bootstrap (one anatomy each; no between-subject variability).

## Anatomies

| anatomy | description | OFC [mm] | breadth x length [mm] | targets | usable cortex [cm2] | OPM dense / matched sites |
|---|---|---|---|---|---|---|
| adult | MNE sample subject (as G2) | 586 | 171 x 211 | 7661 | 1878 | 212 / 96 |
| school-age size (scaled adult) | adult x 85/95 (Jas Table 1 child/adult head radius): 0.8947 | 525 | 153 x 189 | 7507 | 1470 | 171 / 89 |
| 2-year size (scaled adult) | adult x template/adult occipitofrontal circumference: 0.8442 | 495 | 144 x 178 | 7413 | 1290 | 155 / 90 |
| 2-year template | 2-year template (ANTS2-0Years3T), native dimensions | 495 | 140 x 175 | 8134 | 1062 | 151 / 83 |
| 18-month template | 18-month template (ANTS18-0Months3T), native dimensions | 491 | 137 x 172 | 7994 | 975 | 157 / 80 |
| 12-month template | 12-month template (ANTS12-0Months3T), native dimensions | 469 | 133 x 164 | 8088 | 895 | 145 / 83 |

## Placements in the fixed helmet (magnetometer coil centre to scalp)

| anatomy | placement | moved [mm] | min [mm] | median [mm] | feasible |
|---|---|---|---|---|---|
| adult | centred | 0.0 | 23.2 | 29.8 | True |
| adult | top | 5.5 | 20.4 | 28.4 | True |
| adult | back | 5.5 | 20.4 | 28.8 | True |
| adult | x-centred | 5.5 | 20.1 | 28.3 | True |
| adult | top-18mm | 8.0 | 18.0 | 27.8 | True |
| adult | counterfactual | 0.0 | 23.2 | 29.8 | True |
| adult | counterfactual_x-centred | 0.0 | 22.5 | 29.6 | True |
| school-age size (scaled adult) | centred | 0.0 | 32.3 | 41.3 | True |
| school-age size (scaled adult) | top | 22.0 | 20.3 | 34.1 | True |
| school-age size (scaled adult) | back | 16.0 | 20.5 | 38.6 | True |
| school-age size (scaled adult) | x-centred | 22.0 | 20.2 | 34.4 | True |
| school-age size (scaled adult) | top-18mm | 24.5 | 18.0 | 33.6 | True |
| school-age size (scaled adult) | counterfactual | 0.0 | 20.8 | 26.7 | True |
| school-age size (scaled adult) | counterfactual_x-centred | 0.0 | 20.9 | 26.6 | True |
| 2-year size (scaled adult) | centred | 0.0 | 36.6 | 46.9 | True |
| 2-year size (scaled adult) | top | 22.0 | 20.1 | 39.1 | True |
| 2-year size (scaled adult) | back | 21.5 | 20.2 | 44.5 | True |
| 2-year size (scaled adult) | x-centred | 22.0 | 20.0 | 39.4 | True |
| 2-year size (scaled adult) | top-18mm | 24.0 | 18.2 | 38.5 | True |
| 2-year size (scaled adult) | counterfactual | 0.0 | 19.6 | 25.2 | True |
| 2-year size (scaled adult) | counterfactual_x-centred | 0.0 | 19.8 | 25.1 | True |
| 2-year template | centred | 0.0 | 29.0 | 47.3 | True |
| 2-year template | top | 28.0 | 20.2 | 38.5 | True |
| 2-year template | back | 15.0 | 20.0 | 44.3 | True |
| 2-year template | x-centred | 31.0 | 20.3 | 37.5 | True |
| 2-year template | top-18mm | 30.5 | 18.1 | 37.8 | True |
| 2-year template | counterfactual | 0.0 | 18.2 | 34.2 | True |
| 2-year template | counterfactual_x-centred | 0.0 | 18.2 | 27.0 | True |
| 18-month template | centred | 0.0 | 31.8 | 46.0 | True |
| 18-month template | top | 28.0 | 20.4 | 38.9 | True |
| 18-month template | back | 15.5 | 20.2 | 44.8 | True |
| 18-month template | x-centred | 30.5 | 20.1 | 38.3 | True |
| 18-month template | top-18mm | 30.5 | 18.3 | 37.7 | True |
| 18-month template | counterfactual | 0.0 | 18.6 | 31.6 | True |
| 18-month template | counterfactual_x-centred | 0.0 | 18.6 | 27.2 | True |
| 12-month template | centred | 0.0 | 35.5 | 48.9 | True |
| 12-month template | top | 33.5 | 20.2 | 40.4 | True |
| 12-month template | back | 20.0 | 20.3 | 48.1 | True |
| 12-month template | x-centred | 34.5 | 20.1 | 40.2 | True |
| 12-month template | top-18mm | 36.0 | 18.2 | 39.7 | True |
| 12-month template | counterfactual | 0.0 | 18.1 | 29.5 | True |
| 12-month template | counterfactual_x-centred | 0.0 | 18.5 | 27.7 | True |

Link to G2: at the adult's measured (= centred) position the dense/combined detectability ratio is 1.130x as an unweighted median over all targets (G2's headline), 1.143x without the medial wall and +1.25 dB area-weighted without it (the G3B convention).

## D_child, D_adult and Delta (dense OPM; intrinsic + brain noise; detectability dB)

| child anatomy | comparator | D_child | D_adult | Delta | homology |
|---|---|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | +1.49 dB [+1.27, +1.73] | +0.99 dB [+0.75, +1.20] | +0.47 dB [+0.38, +0.58] | vertex |
| school-age size (scaled adult) | Neuromag grad | +2.97 dB [+2.64, +3.29] | +2.21 dB [+2.00, +2.43] | +0.69 dB [+0.55, +0.85] | vertex |
| school-age size (scaled adult) | Neuromag mag | +1.80 dB [+1.55, +2.09] | +1.37 dB [+1.13, +1.61] | +0.44 dB [+0.35, +0.53] | vertex |
| 2-year size (scaled adult) | Neuromag combined | +2.15 dB [+1.84, +2.44] | +0.99 dB [+0.78, +1.20] | +1.13 dB [+0.99, +1.25] | vertex |
| 2-year size (scaled adult) | Neuromag grad | +3.86 dB [+3.50, +4.17] | +2.21 dB [+1.98, +2.43] | +1.53 dB [+1.38, +1.72] | vertex |
| 2-year size (scaled adult) | Neuromag mag | +2.42 dB [+2.08, +2.72] | +1.37 dB [+1.10, +1.61] | +1.03 dB [+0.90, +1.15] | vertex |
| 2-year template | Neuromag combined | +1.85 dB [+1.52, +2.19] | +0.99 dB [+0.76, +1.20] | +0.73 dB [+0.44, +1.16] | parcel |
| 2-year template | Neuromag grad | +3.33 dB [+2.90, +3.75] | +2.21 dB [+1.95, +2.43] | +1.05 dB [+0.69, +1.50] | parcel |
| 2-year template | Neuromag mag | +2.23 dB [+1.91, +2.57] | +1.37 dB [+1.15, +1.62] | +0.82 dB [+0.59, +1.08] | parcel |
| 18-month template | Neuromag combined | +1.90 dB [+1.62, +2.21] | +0.99 dB [+0.77, +1.20] | +0.85 dB [+0.49, +1.27] | parcel |
| 18-month template | Neuromag grad | +3.38 dB [+2.93, +3.85] | +2.21 dB [+2.00, +2.43] | +1.21 dB [+0.82, +1.66] | parcel |
| 18-month template | Neuromag mag | +2.25 dB [+1.94, +2.61] | +1.37 dB [+1.10, +1.61] | +0.67 dB [+0.54, +1.01] | parcel |
| 12-month template | Neuromag combined | +2.06 dB [+1.68, +2.45] | +0.99 dB [+0.76, +1.21] | +1.03 dB [+0.66, +1.44] | parcel |
| 12-month template | Neuromag grad | +3.66 dB [+3.22, +4.25] | +2.21 dB [+1.98, +2.43] | +1.12 dB [+0.92, +2.08] | parcel |
| 12-month template | Neuromag mag | +2.46 dB [+2.01, +2.86] | +1.37 dB [+1.12, +1.61] | +1.01 dB [+0.72, +1.41] | parcel |

Projected condition (room-field subspace removed):

| child anatomy | comparator | D_child | D_adult | Delta |
|---|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | +1.35 dB [+1.07, +1.64] | +0.67 dB [+0.26, +0.96] | +0.90 dB [+0.78, +1.07] |
| school-age size (scaled adult) | Neuromag grad | +2.49 dB [+2.15, +2.83] | +1.51 dB [+1.13, +1.82] | +1.13 dB [+0.98, +1.28] |
| school-age size (scaled adult) | Neuromag mag | +1.86 dB [+1.58, +2.14] | +1.18 dB [+0.78, +1.48] | +0.82 dB [+0.72, +0.96] |
| 2-year size (scaled adult) | Neuromag combined | +2.05 dB [+1.69, +2.39] | +0.67 dB [+0.30, +0.92] | +1.58 dB [+1.48, +1.69] |
| 2-year size (scaled adult) | Neuromag grad | +3.42 dB [+2.94, +3.80] | +1.51 dB [+1.13, +1.81] | +1.97 dB [+1.85, +2.12] |
| 2-year size (scaled adult) | Neuromag mag | +2.50 dB [+2.11, +2.83] | +1.18 dB [+0.81, +1.46] | +1.45 dB [+1.35, +1.57] |
| 2-year template | Neuromag combined | +1.48 dB [+1.04, +1.87] | +0.67 dB [+0.28, +0.94] | +0.74 dB [+0.45, +1.10] |
| 2-year template | Neuromag grad | +2.68 dB [+2.18, +3.22] | +1.51 dB [+1.14, +1.83] | +1.10 dB [+0.77, +1.46] |
| 2-year template | Neuromag mag | +2.02 dB [+1.57, +2.44] | +1.18 dB [+0.80, +1.49] | +0.80 dB [+0.58, +0.91] |
| 18-month template | Neuromag combined | +1.57 dB [+1.19, +1.96] | +0.67 dB [+0.29, +0.94] | +0.91 dB [+0.53, +1.32] |
| 18-month template | Neuromag grad | +2.78 dB [+2.29, +3.28] | +1.51 dB [+1.10, +1.85] | +1.21 dB [+0.83, +1.63] |
| 18-month template | Neuromag mag | +2.06 dB [+1.62, +2.47] | +1.18 dB [+0.81, +1.50] | +0.78 dB [+0.54, +1.02] |
| 12-month template | Neuromag combined | +1.82 dB [+1.32, +2.29] | +0.67 dB [+0.31, +0.94] | +1.17 dB [+0.74, +1.45] |
| 12-month template | Neuromag grad | +3.16 dB [+2.52, +3.68] | +1.51 dB [+1.13, +1.84] | +1.42 dB [+0.96, +2.04] |
| 12-month template | Neuromag mag | +2.37 dB [+1.76, +2.84] | +1.18 dB [+0.81, +1.49] | +1.02 dB [+0.76, +1.49] |

Matched-site OPM (coverage control), intrinsic + brain:

| child anatomy | comparator | D_child | D_adult | Delta |
|---|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | +0.40 dB [+0.23, +0.59] | -0.07 dB [-0.15, +0.03] | +0.46 dB [+0.35, +0.63] |
| school-age size (scaled adult) | Neuromag grad | +1.65 dB [+1.45, +1.92] | +0.85 dB [+0.67, +0.99] | +0.67 dB [+0.50, +0.88] |
| school-age size (scaled adult) | Neuromag mag | +0.68 dB [+0.49, +0.83] | +0.21 dB [+0.06, +0.36] | +0.41 dB [+0.31, +0.56] |
| 2-year size (scaled adult) | Neuromag combined | +1.13 dB [+0.86, +1.37] | -0.07 dB [-0.17, +0.04] | +1.22 dB [+1.10, +1.36] |
| 2-year size (scaled adult) | Neuromag grad | +2.60 dB [+2.28, +2.89] | +0.85 dB [+0.69, +0.99] | +1.62 dB [+1.43, +1.81] |
| 2-year size (scaled adult) | Neuromag mag | +1.40 dB [+1.14, +1.65] | +0.21 dB [+0.06, +0.35] | +1.13 dB [+0.99, +1.26] |
| 2-year template | Neuromag combined | +0.55 dB [+0.30, +0.82] | -0.07 dB [-0.16, +0.02] | +0.90 dB [+0.48, +1.06] |
| 2-year template | Neuromag grad | +1.95 dB [+1.57, +2.39] | +0.85 dB [+0.69, +1.00] | +1.21 dB [+0.75, +1.63] |
| 2-year template | Neuromag mag | +0.88 dB [+0.58, +1.16] | +0.21 dB [+0.05, +0.37] | +0.79 dB [+0.48, +0.98] |
| 18-month template | Neuromag combined | +0.57 dB [+0.33, +0.83] | -0.07 dB [-0.18, +0.02] | +0.85 dB [+0.46, +1.15] |
| 18-month template | Neuromag grad | +1.99 dB [+1.63, +2.35] | +0.85 dB [+0.69, +1.00] | +1.05 dB [+0.79, +1.64] |
| 18-month template | Neuromag mag | +0.87 dB [+0.61, +1.08] | +0.21 dB [+0.05, +0.38] | +0.77 dB [+0.32, +0.90] |
| 12-month template | Neuromag combined | +0.74 dB [+0.47, +1.03] | -0.07 dB [-0.18, +0.03] | +1.07 dB [+0.54, +1.39] |
| 12-month template | Neuromag grad | +2.33 dB [+1.88, +2.86] | +0.85 dB [+0.68, +1.02] | +1.35 dB [+0.93, +2.34] |
| 12-month template | Neuromag mag | +1.07 dB [+0.75, +1.42] | +0.21 dB [+0.06, +0.35] | +0.87 dB [+0.46, +1.34] |

## What drives Delta: placement and helmet fit (dense OPM vs Neuromag combined, intrinsic + brain)

Each child placement is compared with the adult at the same rule (the adult's counterfactual helmet has factor 1).

| child anatomy | top (primary) | centred | x-centred | top-18mm | back | counterfactual | counterfactual, x-centred |
|---|---|---|---|---|---|---|---|
| school-age size (scaled adult) | +0.47 dB [+0.38, +0.58] | +1.15 dB [+1.03, +1.24] | +0.47 dB [+0.38, +0.57] | +0.46 dB [+0.37, +0.57] | +0.89 dB [+0.75, +1.05] | -0.17 dB [-0.21, -0.13] | -0.17 dB [-0.22, -0.13] |
| 2-year size (scaled adult) | +1.13 dB [+0.99, +1.25] | +1.73 dB [+1.57, +1.86] | +1.12 dB [+0.99, +1.24] | +1.16 dB [+1.02, +1.27] | +1.30 dB [+0.99, +1.53] | -0.30 dB [-0.35, -0.25] | -0.30 dB [-0.36, -0.24] |
| 2-year template | +0.73 dB [+0.44, +1.16] | +1.79 dB [+1.47, +2.30] | +0.68 dB [+0.39, +0.79] | +0.70 dB [+0.45, +1.01] | +1.39 dB [+1.24, +1.63] | +0.28 dB [+0.05, +0.65] | -0.37 dB [-0.51, +0.02] |
| 18-month template | +0.85 dB [+0.49, +1.27] | +1.93 dB [+1.74, +2.12] | +0.69 dB [+0.47, +0.97] | +0.83 dB [+0.49, +1.15] | +1.63 dB [+1.24, +1.80] | +0.20 dB [-0.05, +0.58] | -0.29 dB [-0.45, +0.01] |
| 12-month template | +1.03 dB [+0.66, +1.44] | +2.49 dB [+2.14, +2.91] | +0.99 dB [+0.66, +1.38] | +0.95 dB [+0.54, +1.40] | +1.68 dB [+1.19, +2.18] | +0.08 dB [-0.23, +0.22] | -0.12 dB [-0.41, +0.10] |

Counterfactual Delta by comparator (helmet scaled with the head; the dependence on the comparator points to the SQUID side of the change):

| child anatomy | comparator | counterfactual | counterfactual, x-centred |
|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | -0.17 dB [-0.21, -0.13] | -0.17 dB [-0.22, -0.13] |
| school-age size (scaled adult) | Neuromag grad | -0.32 dB [-0.36, -0.27] | -0.32 dB [-0.38, -0.26] |
| school-age size (scaled adult) | Neuromag mag | -0.12 dB [-0.16, -0.08] | -0.12 dB [-0.16, -0.08] |
| 2-year size (scaled adult) | Neuromag combined | -0.30 dB [-0.35, -0.25] | -0.30 dB [-0.36, -0.24] |
| 2-year size (scaled adult) | Neuromag grad | -0.52 dB [-0.59, -0.45] | -0.52 dB [-0.58, -0.45] |
| 2-year size (scaled adult) | Neuromag mag | -0.21 dB [-0.26, -0.16] | -0.21 dB [-0.26, -0.17] |
| 2-year template | Neuromag combined | +0.28 dB [+0.05, +0.65] | -0.37 dB [-0.51, +0.02] |
| 2-year template | Neuromag grad | +0.22 dB [+0.03, +0.50] | -0.50 dB [-0.67, -0.31] |
| 2-year template | Neuromag mag | +0.38 dB [+0.27, +0.61] | -0.16 dB [-0.31, +0.15] |
| 18-month template | Neuromag combined | +0.20 dB [-0.05, +0.58] | -0.29 dB [-0.45, +0.01] |
| 18-month template | Neuromag grad | +0.25 dB [+0.12, +0.48] | -0.40 dB [-0.57, -0.04] |
| 18-month template | Neuromag mag | +0.28 dB [+0.10, +0.45] | -0.23 dB [-0.31, +0.11] |
| 12-month template | Neuromag combined | +0.08 dB [-0.23, +0.22] | -0.12 dB [-0.41, +0.10] |
| 12-month template | Neuromag grad | -0.13 dB [-0.34, +0.01] | -0.32 dB [-0.60, -0.07] |
| 12-month template | Neuromag mag | +0.09 dB [-0.01, +0.31] | -0.00 dB [-0.15, +0.16] |

## Absolute detectability (median 20 log10 d of a 10-nAm dipole, intrinsic + brain, primary placement)

| anatomy | OPM dense | OPM matched | Neuromag combined | Neuromag grad | Neuromag mag |
|---|---|---|---|---|---|
| adult | -1.60 | -2.80 | -2.66 | -3.63 | -3.00 |
| school-age size (scaled adult) | -0.66 | -1.74 | -2.15 | -3.34 | -2.42 |
| 2-year size (scaled adult) | -0.19 | -1.18 | -2.27 | -3.66 | -2.52 |
| 2-year template | +0.93 | -0.40 | -1.19 | -2.70 | -1.52 |
| 18-month template | +1.08 | -0.29 | -1.02 | -2.66 | -1.34 |
| 12-month template | +1.87 | +0.46 | -0.62 | -2.23 | -0.96 |

Vertex-wise change from the adult (scaled controls; same vertex): both systems gain, the OPM more.

| child anatomy | OPM dense | Neuromag combined | Neuromag grad | Neuromag mag |
|---|---|---|---|---|
| school-age size (scaled adult) | +1.15 | +0.65 | +0.45 | +0.69 |
| 2-year size (scaled adult) | +1.74 | +0.54 | +0.10 | +0.64 |

## Channel count: the adult's dense array subsampled to each child's site count

| child anatomy | sites | D_child | D_adult, subsampled | Delta at equal channel count |
|---|---|---|---|---|
| school-age size (scaled adult) | 171 | +1.49 dB [+1.31, +1.71] | +0.66 dB [+0.49, +0.85] | +0.72 dB [+0.58, +0.84] |
| 2-year size (scaled adult) | 155 | +2.15 dB [+1.82, +2.49] | +0.52 dB [+0.35, +0.64] | +1.50 dB [+1.33, +1.66] |
| 2-year template | 151 | +1.85 dB [+1.56, +2.19] | +0.49 dB [+0.34, +0.66] | +1.36 dB [+0.90, +1.69] |
| 18-month template | 157 | +1.90 dB [+1.65, +2.20] | +0.54 dB [+0.40, +0.68] | +1.37 dB [+0.88, +1.61] |
| 12-month template | 145 | +2.06 dB [+1.68, +2.43] | +0.41 dB [+0.26, +0.56] | +1.49 dB [+1.25, +2.32] |

## Delta by depth stratum (dense OPM vs Neuromag combined, intrinsic + brain)

| child anatomy | depth [mm] | n child / adult | D_child | D_adult | Delta [95 % CI] |
|---|---|---|---|---|---|
| school-age size (scaled adult) | 0-10 | 0 / 0 | sparse | | |
| school-age size (scaled adult) | 10-15 | 450 / 222 | +3.85 | +3.72 | +0.12 [-0.39, +0.60] |
| school-age size (scaled adult) | 15-20 | 1604 / 1221 | +2.82 | +2.42 | +0.39 [+0.16, +0.63] |
| school-age size (scaled adult) | 20-25 | 1604 / 1572 | +1.76 | +1.54 | +0.23 [+0.05, +0.41] |
| school-age size (scaled adult) | 25-30 | 1144 / 1289 | +1.03 | +0.88 | +0.16 [+0.02, +0.30] |
| school-age size (scaled adult) | 30-40 | 1399 / 1539 | +0.59 | +0.40 | +0.19 [+0.07, +0.29] |
| school-age size (scaled adult) | 40-50 | 633 / 913 | +0.43 | +0.17 | +0.25 [+0.13, +0.38] |
| school-age size (scaled adult) | 50-60 | 115 / 321 | +0.78 | +0.08 | +0.70 [-0.32, +0.87] |
| school-age size (scaled adult) | 60-90 | 0 / 26 | sparse | | |
| 2-year size (scaled adult) | 0-10 | 2 / 0 | sparse | | |
| 2-year size (scaled adult) | 10-15 | 610 / 222 | +4.94 | +3.72 | +1.21 [+0.74, +1.67] |
| 2-year size (scaled adult) | 15-20 | 1768 / 1221 | +3.49 | +2.42 | +1.06 [+0.83, +1.32] |
| 2-year size (scaled adult) | 20-25 | 1610 / 1572 | +2.19 | +1.54 | +0.66 [+0.47, +0.85] |
| 2-year size (scaled adult) | 25-30 | 1055 / 1289 | +1.28 | +0.88 | +0.41 [+0.22, +0.58] |
| 2-year size (scaled adult) | 30-40 | 1301 / 1539 | +0.78 | +0.40 | +0.38 [+0.25, +0.53] |
| 2-year size (scaled adult) | 40-50 | 463 / 913 | +0.63 | +0.17 | +0.46 [+0.32, +0.65] |
| 2-year size (scaled adult) | 50-60 | 46 / 321 | +1.13 | +0.08 | +1.05 [-0.14, +1.22] |
| 2-year size (scaled adult) | 60-90 | 0 / 26 | sparse | | |
| 2-year template | 0-10 | 26 / 0 | sparse | | |
| 2-year template | 10-15 | 1488 / 222 | +4.08 | +3.72 | +0.35 [-0.15, +0.95] |
| 2-year template | 15-20 | 1666 / 1221 | +2.63 | +2.42 | +0.21 [-0.15, +0.58] |
| 2-year template | 20-25 | 1410 / 1572 | +1.65 | +1.54 | +0.12 [-0.10, +0.35] |
| 2-year template | 25-30 | 980 / 1289 | +1.04 | +0.88 | +0.16 [-0.02, +0.35] |
| 2-year template | 30-40 | 1175 / 1539 | +0.73 | +0.40 | +0.33 [+0.17, +0.60] |
| 2-year template | 40-50 | 502 / 913 | +0.55 | +0.17 | +0.37 [+0.26, +0.74] |
| 2-year template | 50-60 | 229 / 321 | +1.20 | +0.08 | +1.12 [+0.31, +1.42] |
| 2-year template | 60-90 | 83 / 26 | +2.23 | +0.00 | +2.23 [+1.84, +2.43] |
| 18-month template | 0-10 | 32 / 0 | sparse | | |
| 18-month template | 10-15 | 744 / 222 | +4.21 | +3.72 | +0.49 [-0.09, +1.18] |
| 18-month template | 15-20 | 1788 / 1221 | +2.95 | +2.42 | +0.52 [+0.12, +1.08] |
| 18-month template | 20-25 | 1539 / 1572 | +1.92 | +1.54 | +0.39 [+0.07, +0.80] |
| 18-month template | 25-30 | 1161 / 1289 | +1.31 | +0.88 | +0.44 [+0.20, +0.66] |
| 18-month template | 30-40 | 1290 / 1539 | +0.85 | +0.40 | +0.45 [+0.28, +0.68] |
| 18-month template | 40-50 | 523 / 913 | +0.54 | +0.17 | +0.37 [+0.27, +0.69] |
| 18-month template | 50-60 | 221 / 321 | +1.26 | +0.08 | +1.18 [+0.40, +1.44] |
| 18-month template | 60-90 | 65 / 26 | +2.15 | +0.00 | +2.15 [+1.74, +2.22] |
| 12-month template | 0-10 | 78 / 0 | sparse | | |
| 12-month template | 10-15 | 1744 / 222 | +4.55 | +3.72 | +0.83 [+0.12, +1.55] |
| 12-month template | 15-20 | 1768 / 1221 | +2.90 | +2.42 | +0.48 [-0.08, +0.90] |
| 12-month template | 20-25 | 1452 / 1572 | +1.78 | +1.54 | +0.24 [-0.09, +0.56] |
| 12-month template | 25-30 | 928 / 1289 | +1.09 | +0.88 | +0.21 [+0.04, +0.43] |
| 12-month template | 30-40 | 1029 / 1539 | +0.71 | +0.40 | +0.31 [+0.15, +0.51] |
| 12-month template | 40-50 | 390 / 913 | +0.61 | +0.17 | +0.44 [+0.29, +0.61] |
| 12-month template | 50-60 | 184 / 321 | +1.15 | +0.08 | +1.07 [+0.56, +1.30] |
| 12-month template | 60-90 | 19 / 26 | +1.53 | +0.00 | +1.53 [+1.11, +1.56] |

Scaled controls, vertex-wise (homologous) Delta by the adult's depth:

| child anatomy | adult depth [mm] | n | Delta [95 % CI] |
|---|---|---|---|
| school-age size (scaled adult) | 0-10 | 0 | sparse |
| school-age size (scaled adult) | 10-15 | 175 | +0.73 [+0.24, +1.16] |
| school-age size (scaled adult) | 15-20 | 1134 | +0.79 [+0.60, +1.04] |
| school-age size (scaled adult) | 20-25 | 1564 | +0.70 [+0.57, +0.84] |
| school-age size (scaled adult) | 25-30 | 1284 | +0.52 [+0.43, +0.62] |
| school-age size (scaled adult) | 30-40 | 1533 | +0.31 [+0.24, +0.38] |
| school-age size (scaled adult) | 40-50 | 912 | +0.22 [+0.13, +0.35] |
| school-age size (scaled adult) | 50-60 | 321 | +0.35 [+0.15, +0.54] |
| school-age size (scaled adult) | 60-90 | 26 | +0.91 [-0.58, +1.28] |
| 2-year size (scaled adult) | 0-10 | 0 | sparse |
| 2-year size (scaled adult) | 10-15 | 152 | +2.25 [+1.67, +2.53] |
| 2-year size (scaled adult) | 15-20 | 1082 | +1.89 [+1.76, +2.03] |
| 2-year size (scaled adult) | 20-25 | 1551 | +1.56 [+1.43, +1.68] |
| 2-year size (scaled adult) | 25-30 | 1283 | +1.15 [+1.04, +1.26] |
| 2-year size (scaled adult) | 30-40 | 1532 | +0.68 [+0.61, +0.78] |
| 2-year size (scaled adult) | 40-50 | 909 | +0.45 [+0.35, +0.58] |
| 2-year size (scaled adult) | 50-60 | 320 | +0.49 [+0.30, +0.70] |
| 2-year size (scaled adult) | 60-90 | 26 | +1.10 [-0.26, +1.62] |

## Delta by orientation stratum (0 deg = radial to the inner skull; dense OPM vs Neuromag combined)

| child anatomy | orientation [deg] | n child / adult | D_child | D_adult | Delta [95 % CI] |
|---|---|---|---|---|---|
| school-age size (scaled adult) | 0-30 | 676 / 707 | +1.27 | +0.57 | +0.70 [+0.46, +0.98] |
| school-age size (scaled adult) | 30-60 | 2297 / 2346 | +1.27 | +0.76 | +0.51 [+0.21, +0.80] |
| school-age size (scaled adult) | 60-90.1 | 3976 / 4050 | +1.66 | +1.19 | +0.47 [+0.13, +0.83] |
| 2-year size (scaled adult) | 0-30 | 665 / 707 | +1.80 | +0.57 | +1.23 [+0.93, +1.47] |
| 2-year size (scaled adult) | 30-60 | 2269 / 2346 | +1.86 | +0.76 | +1.09 [+0.72, +1.49] |
| 2-year size (scaled adult) | 60-90.1 | 3921 / 4050 | +2.43 | +1.19 | +1.24 [+0.80, +1.63] |
| 2-year template | 0-30 | 1219 / 707 | +2.78 | +0.57 | +2.22 [+1.58, +2.82] |
| 2-year template | 30-60 | 2293 / 2346 | +1.76 | +0.76 | +1.00 [+0.62, +1.37] |
| 2-year template | 60-90.1 | 4047 / 4050 | +1.72 | +1.19 | +0.53 [+0.18, +0.90] |
| 18-month template | 0-30 | 1026 / 707 | +2.57 | +0.57 | +2.00 [+1.43, +2.48] |
| 18-month template | 30-60 | 2398 / 2346 | +1.91 | +0.76 | +1.15 [+0.73, +1.55] |
| 18-month template | 60-90.1 | 3939 / 4050 | +1.78 | +1.19 | +0.59 [+0.27, +0.99] |
| 12-month template | 0-30 | 1294 / 707 | +3.25 | +0.57 | +2.68 [+1.77, +3.65] |
| 12-month template | 30-60 | 2461 / 2346 | +1.90 | +0.76 | +1.14 [+0.69, +1.63] |
| 12-month template | 60-90.1 | 3837 / 4050 | +1.89 | +1.19 | +0.70 [+0.31, +1.13] |

## Placement, counterfactual helmet and sensitivity (median D, dense OPM vs Neuromag combined, intrinsic + brain)

| anatomy | centred | top | back | x+5mm | x-5mm | y+5mm | y-5mm | pitch+10deg | pitch-10deg | roll+5deg | roll-5deg | x-centred | top-18mm | counterfactual | counterfactual_x-centred |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| adult | +1.25 | +0.99 | +1.11 | +1.12 | +1.26 | +1.23 | +0.97 | +0.96 | +1.41 | +1.25 | +1.15 | +1.00 | +0.86 | +1.25 | +1.25 |
| school-age size (scaled adult) | +2.47 | +1.49 | +1.99 | +1.47 | +1.51 | +1.55 | +1.45 | +1.48 | +2.52 | +1.45 | +1.69 | +1.49 | +1.35 | +0.93 | +0.94 |
| 2-year size (scaled adult) | +3.06 | +2.15 | +2.36 | +2.01 | +2.16 | +2.24 | +2.07 | +1.70 | +3.06 | +1.94 | +2.53 | +2.16 | +2.05 | +0.75 | +0.75 |
| 2-year template | +3.08 | +1.85 | +2.51 | +2.07 | +1.75 | +1.84 | +1.85 | +1.78 | +2.01 | +1.79 | +2.05 | +1.76 | +1.70 | +1.56 | +1.01 |
| 18-month template | +3.13 | +1.90 | +2.49 | +2.05 | +1.81 | +1.87 | +1.82 | +1.81 | +2.04 | +1.81 | +2.04 | +1.81 | +1.75 | +1.47 | +1.03 |
| 12-month template | +3.66 | +2.06 | +2.80 | +2.14 | +2.04 | +2.09 | +2.08 | +1.99 | +2.23 | +2.02 | +2.12 | +1.99 | +1.86 | +1.18 | +1.01 |

| anatomy | opm_asd_7fT | opm_asd_10fT | opm_asd_15fT | opm_asd_20fT | opm_asd_30fT | background_x0.5 | background_x2 | bem1 |
|---|---|---|---|---|---|---|---|---|
| adult | +1.99 | +1.53 | +0.99 | +0.57 | -0.04 | +0.97 | +0.99 | +1.02 |
| school-age size (scaled adult) | +2.48 | +2.08 | +1.49 | +1.02 | +0.34 | +1.48 | +1.47 | +1.45 |
| 2-year size (scaled adult) | +3.16 | +2.73 | +2.15 | +1.69 | +0.99 | +2.17 | +2.09 | +2.12 |
| 2-year template | +2.82 | +2.40 | +1.85 | +1.38 | +0.61 | +1.86 | +1.79 | +1.80 |
| 18-month template | +2.98 | +2.54 | +1.90 | +1.40 | +0.61 | +1.86 | +1.91 | +1.81 |
| 12-month template | +3.09 | +2.66 | +2.06 | +1.55 | +0.77 | +2.07 | +2.00 | +2.02 |

Difference of these medians, child minus adult, with the same variant applied to both (a sensitivity of the medians, not the paired Delta estimator; the background variants scale the adult too):

| child anatomy | opm_asd_7fT | opm_asd_10fT | opm_asd_15fT | opm_asd_20fT | opm_asd_30fT | background_x0.5 | background_x2 | bem1 |
|---|---|---|---|---|---|---|---|---|
| school-age size (scaled adult) | +0.49 | +0.55 | +0.49 | +0.45 | +0.38 | +0.51 | +0.48 | +0.43 |
| 2-year size (scaled adult) | +1.17 | +1.20 | +1.16 | +1.12 | +1.03 | +1.20 | +1.10 | +1.10 |
| 2-year template | +0.82 | +0.87 | +0.85 | +0.81 | +0.64 | +0.89 | +0.80 | +0.78 |
| 18-month template | +0.99 | +1.01 | +0.91 | +0.83 | +0.65 | +0.89 | +0.91 | +0.79 |
| 12-month template | +1.09 | +1.13 | +1.07 | +0.98 | +0.81 | +1.10 | +1.01 | +1.00 |

## Regions that gain or lose with the placement (median D by lobe, dense OPM vs Neuromag combined, intrinsic + brain)

| anatomy | placement | frontal | parietal | temporal | occipital | cingulate | insula |
|---|---|---|---|---|---|---|---|
| adult | centred | +1.85 | +1.13 | +1.17 | +0.74 | +0.27 | +0.58 |
| adult | top | +1.49 | +0.84 | +1.02 | +0.63 | +0.15 | +0.37 |
| adult | x-centred | +1.52 | +0.87 | +1.05 | +0.63 | +0.16 | +0.37 |
| adult | back | +1.99 | +0.94 | +1.04 | +0.39 | +0.26 | +0.47 |
| adult | counterfactual | +1.85 | +1.13 | +1.17 | +0.74 | +0.27 | +0.58 |
| school-age size (scaled adult) | centred | +3.09 | +2.55 | +2.51 | +1.99 | +0.89 | +1.33 |
| school-age size (scaled adult) | top | +1.72 | +1.33 | +1.88 | +1.50 | +0.37 | +0.80 |
| school-age size (scaled adult) | x-centred | +1.73 | +1.34 | +1.89 | +1.49 | +0.38 | +0.79 |
| school-age size (scaled adult) | back | +3.35 | +1.87 | +1.97 | +0.74 | +0.77 | +1.12 |
| school-age size (scaled adult) | counterfactual | +1.35 | +0.92 | +0.90 | +0.61 | +0.31 | +0.50 |
| 2-year size (scaled adult) | centred | +3.78 | +3.05 | +3.20 | +2.37 | +1.12 | +1.70 |
| 2-year size (scaled adult) | top | +2.58 | +2.03 | +2.54 | +1.93 | +0.59 | +1.23 |
| 2-year size (scaled adult) | x-centred | +2.59 | +2.04 | +2.53 | +1.93 | +0.59 | +1.23 |
| 2-year size (scaled adult) | back | +4.04 | +2.12 | +2.38 | +0.77 | +0.92 | +1.39 |
| 2-year size (scaled adult) | counterfactual | +1.19 | +0.68 | +0.73 | +0.39 | +0.21 | +0.43 |
| 2-year template | centred | +3.71 | +3.02 | +3.29 | +2.75 | +1.10 | +1.97 |
| 2-year template | top | +2.02 | +1.39 | +2.63 | +2.08 | +0.61 | +1.36 |
| 2-year template | x-centred | +1.95 | +1.30 | +2.64 | +2.02 | +0.57 | +1.29 |
| 2-year template | back | +3.89 | +2.20 | +2.77 | +1.25 | +1.07 | +1.75 |
| 2-year template | counterfactual | +2.16 | +1.25 | +1.72 | +1.21 | +0.67 | +1.17 |
| 18-month template | centred | +3.86 | +2.80 | +3.51 | +2.59 | +1.08 | +1.88 |
| 18-month template | top | +2.25 | +1.34 | +2.90 | +1.92 | +0.56 | +1.39 |
| 18-month template | x-centred | +2.15 | +1.24 | +2.86 | +1.86 | +0.54 | +1.28 |
| 18-month template | back | +4.03 | +1.99 | +2.82 | +1.20 | +1.05 | +1.78 |
| 18-month template | counterfactual | +2.14 | +1.03 | +1.76 | +1.01 | +0.58 | +1.09 |
| 12-month template | centred | +4.51 | +3.59 | +3.83 | +3.28 | +1.23 | +1.84 |
| 12-month template | top | +2.46 | +1.44 | +3.06 | +2.34 | +0.62 | +1.26 |
| 12-month template | x-centred | +2.42 | +1.37 | +3.02 | +2.31 | +0.60 | +1.24 |
| 12-month template | back | +4.69 | +2.58 | +2.93 | +1.11 | +1.06 | +1.64 |
| 12-month template | counterfactual | +2.05 | +0.85 | +1.26 | +0.70 | +0.48 | +0.72 |

## Source-to-sensor distance (median, mm) by depth below the scalp

| anatomy | sensors | 0-10 | 10-15 | 15-20 | 20-25 | 25-30 | 30-40 | 40-50 | 50-60 | 60-90 |
|---|---|---|---|---|---|---|---|---|---|---|
| adult | squid:top | - | 45 | 48 | 52 | 57 | 64 | 74 | 82 | 87 |
| adult | opm_dense | - | 22 | 26 | 31 | 36 | 43 | 53 | 63 | 69 |
| school-age size (scaled adult) | squid:top | - | 48 | 52 | 56 | 61 | 69 | 75 | 78 | - |
| school-age size (scaled adult) | opm_dense | - | 22 | 26 | 31 | 35 | 43 | 52 | 61 | - |
| 2-year size (scaled adult) | squid:top | 57 | 54 | 57 | 62 | 67 | 75 | 80 | 84 | - |
| 2-year size (scaled adult) | opm_dense | 17 | 22 | 26 | 31 | 36 | 43 | 53 | 59 | - |
| 2-year template | squid:top | 45 | 47 | 50 | 56 | 60 | 69 | 74 | 83 | 89 |
| 2-year template | opm_dense | 19 | 21 | 26 | 30 | 35 | 42 | 52 | 61 | 70 |
| 18-month template | squid:top | 44 | 48 | 51 | 56 | 61 | 70 | 75 | 83 | 90 |
| 18-month template | opm_dense | 18 | 22 | 26 | 30 | 35 | 42 | 52 | 62 | 69 |
| 12-month template | squid:top | 50 | 47 | 51 | 57 | 60 | 70 | 73 | 83 | 85 |
| 12-month template | opm_dense | 19 | 21 | 26 | 30 | 35 | 42 | 51 | 62 | 68 |

## Noise composition (median per-channel RMS in the band; magnetometers and OPM in fT, gradiometers in fT/cm)

| anatomy | channels | brain background | intrinsic |
|---|---|---|---|
| adult | squid:top/grad | 41.2 | 21.3 |
| adult | squid:top/mag | 202.6 | 20.7 |
| adult | squid:centred/grad | 37.1 | 21.3 |
| adult | squid:centred/mag | 192.1 | 20.7 |
| adult | opm_dense/mag | 498.4 | 88.8 |
| adult | opm_matched/mag | 516.0 | 88.8 |
| school-age size (scaled adult) | squid:top/grad | 29.8 | 21.3 |
| school-age size (scaled adult) | squid:top/mag | 149.8 | 20.7 |
| school-age size (scaled adult) | squid:centred/grad | 21.4 | 21.3 |
| school-age size (scaled adult) | squid:centred/mag | 131.0 | 20.7 |
| school-age size (scaled adult) | opm_dense/mag | 521.2 | 88.8 |
| school-age size (scaled adult) | opm_matched/mag | 530.7 | 88.8 |
| 2-year size (scaled adult) | squid:top/grad | 23.4 | 21.3 |
| 2-year size (scaled adult) | squid:top/mag | 124.9 | 20.7 |
| 2-year size (scaled adult) | squid:centred/grad | 16.6 | 21.3 |
| 2-year size (scaled adult) | squid:centred/mag | 108.0 | 20.7 |
| 2-year size (scaled adult) | opm_dense/mag | 541.5 | 88.8 |
| 2-year size (scaled adult) | opm_matched/mag | 544.3 | 88.8 |
| 2-year template | squid:top/grad | 23.4 | 21.3 |
| 2-year template | squid:top/mag | 119.6 | 20.7 |
| 2-year template | squid:centred/grad | 16.4 | 21.3 |
| 2-year template | squid:centred/mag | 100.1 | 20.7 |
| 2-year template | opm_dense/mag | 488.4 | 88.8 |
| 2-year template | opm_matched/mag | 483.8 | 88.8 |
| 18-month template | squid:top/grad | 21.1 | 21.3 |
| 18-month template | squid:top/mag | 110.0 | 20.7 |
| 18-month template | squid:centred/grad | 15.7 | 21.3 |
| 18-month template | squid:centred/mag | 94.2 | 20.7 |
| 18-month template | opm_dense/mag | 447.5 | 88.8 |
| 18-month template | opm_matched/mag | 462.6 | 88.8 |
| 12-month template | squid:top/grad | 18.5 | 21.3 |
| 12-month template | squid:top/mag | 98.2 | 20.7 |
| 12-month template | squid:centred/grad | 14.1 | 21.3 |
| 12-month template | squid:centred/mag | 88.1 | 20.7 |
| 12-month template | opm_dense/mag | 473.7 | 88.8 |
| 12-month template | opm_matched/mag | 469.6 | 88.8 |

## Usefulness (dense OPM vs Neuromag combined, intrinsic + brain): share of usable cortical area

A source counts as usable when its detectability reaches 5 at the reference moment (an operational choice, not a clinical standard).

| anatomy | moment [nAm] | both | OPM only | SQUID only | neither |
|---|---|---|---|---|---|
| adult | 20 | 0.00 | 0.04 | 0.00 | 0.96 |
| adult | 50 | 0.35 | 0.07 | 0.00 | 0.57 |
| adult | 100 | 0.66 | 0.02 | 0.00 | 0.32 |
| adult | 200 | 0.89 | 0.01 | 0.00 | 0.10 |
| school-age size (scaled adult) | 20 | 0.00 | 0.06 | 0.00 | 0.94 |
| school-age size (scaled adult) | 50 | 0.37 | 0.10 | 0.00 | 0.53 |
| school-age size (scaled adult) | 100 | 0.69 | 0.04 | 0.00 | 0.28 |
| school-age size (scaled adult) | 200 | 0.91 | 0.02 | 0.00 | 0.07 |
| 2-year size (scaled adult) | 20 | 0.00 | 0.08 | 0.00 | 0.92 |
| 2-year size (scaled adult) | 50 | 0.35 | 0.14 | 0.00 | 0.51 |
| 2-year size (scaled adult) | 100 | 0.69 | 0.05 | 0.00 | 0.26 |
| 2-year size (scaled adult) | 200 | 0.91 | 0.02 | 0.00 | 0.06 |
| 2-year template | 20 | 0.01 | 0.11 | 0.00 | 0.88 |
| 2-year template | 50 | 0.42 | 0.12 | 0.00 | 0.46 |
| 2-year template | 100 | 0.75 | 0.05 | 0.00 | 0.20 |
| 2-year template | 200 | 0.92 | 0.02 | 0.00 | 0.06 |
| 18-month template | 20 | 0.01 | 0.12 | 0.00 | 0.87 |
| 18-month template | 50 | 0.43 | 0.12 | 0.00 | 0.45 |
| 18-month template | 100 | 0.76 | 0.05 | 0.00 | 0.19 |
| 18-month template | 200 | 0.93 | 0.02 | 0.00 | 0.05 |
| 12-month template | 20 | 0.02 | 0.15 | 0.00 | 0.83 |
| 12-month template | 50 | 0.46 | 0.13 | 0.00 | 0.41 |
| 12-month template | 100 | 0.79 | 0.05 | 0.00 | 0.16 |
| 12-month template | 200 | 0.94 | 0.02 | 0.00 | 0.05 |

## Notes

- D = 20 log10(d_OPM / d_SQUID) of a 10-nAm cortical-normal dipole (known-topography detectability with the oracle noise covariance; independent of the moment). Delta = D_child - D_adult. A positive Delta is an increase in relative OPM performance under these matching assumptions; it does not by itself mean that OPM beats SQUID in the child.
- Scaled controls: the adult's vertices, so Delta is vertex-wise. The absolute 4-mm usable-source rule drops 154 and 248 superficial adult targets in the scaled copies, so D_child and D_adult are medians over slightly different target sets while Delta uses the common vertices. The templates: no vertex correspondence; Delta is computed per Desikan-Killiany parcel and per declared depth/orientation stratum from area-weighted medians.
- Intervals: bootstrap over parcels of one anatomy (or of each anatomy, for between-anatomy strata); they do not include between-subject variability. 3 average templates of one database (2-year template, 18-month template, 12-month template) are not a population: template results are conditional simulations.
- Every child array uses the adult's conventions: background moment variance per unit cortical area, room field, intrinsic noise, sensor sizes and the 3-layer BEM conductivities; only geometry changes. Both systems' detectability rises in the smaller heads, the OPM's more (absolute detectability table), by different routes: the on-scalp OPM sees more signal from a cortex that is closer in absolute terms at about the same brain noise, while the SQUIDs' brain noise falls (the cortex is farther from the fixed helmet and, with the background fixed per unit area, smaller) more than their signal. The templates' averaged white surfaces are smoother than an individual cortex (usable area 1,062, 975, 895 cm^2 for the 2-year template, 18-month template, 12-month template vs 1,878 cm^2 for the adult), which lowers their background power and their patch cancellation further; scaling the background variance x0.5 or x2 leaves D_child almost unchanged.
- Placements are chosen from the scalp and helmet geometry only. Under the adult's measured pose a head with other fiducials need not be centred laterally; 'x-centred' shifts each head along device x to equal left/right median gaps before the top contact (shift: adult -1.5 mm, school-age size (scaled adult) -0.5 mm, 2-year size (scaled adult) -0.5 mm, 2-year template -6.5 mm, 18-month template -5.5 mm, 12-month template -2.5 mm; negative = to the left), and 'counterfactual_x-centred' scales the helmet about that laterally centred head. The counterfactual helmet (scaled with the head) is a mechanistic control, not a pediatric SQUID system.
- Targets on the medial wall (FreeSurfer 'unknown': the cut through the corpus callosum and midbrain, not cortex) are left out of every summary; they would otherwise dominate the deepest strata.
