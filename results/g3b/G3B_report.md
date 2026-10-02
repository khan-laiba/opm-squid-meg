# G3B: fixed adult Neuromag helmet versus head-adaptive OPM on smaller heads (NEW)

Code commit: see `g3b_summary.json` (provenance). Configuration: `configs/g3b_pediatric.toml`. Primary placement: `top` (head raised to 20-mm contact). Metric: known-topography detectability in dB; D = OPM - SQUID; Delta = D_child - D_adult. Intervals: parcel bootstrap (one anatomy each; no between-subject variability).

## Anatomies

| anatomy | description | OFC [mm] | breadth x length [mm] | targets | usable cortex [cm2] | OPM dense / matched sites |
|---|---|---|---|---|---|---|
| adult | MNE sample subject (as G2) | 586 | 171 x 211 | 7661 | 1878 | 205 / 94 |
| school-age size (scaled adult) | adult x 85/95 (Jas Table 1 child/adult head radius): 0.8947 | 525 | 153 x 189 | 7507 | 1470 | 172 / 89 |
| 2-year size (scaled adult) | adult x template/adult occipitofrontal circumference: 0.8442 | 495 | 144 x 178 | 7413 | 1290 | 155 / 90 |
| 2-year template | 2-year template (ANTS2-0Years3T), native dimensions | 495 | 140 x 175 | 8134 | 1062 | 151 / 83 |
| 18-month template | 18-month template (ANTS18-0Months3T), native dimensions | 491 | 137 x 172 | 7994 | 975 | 157 / 80 |
| 12-month template | 12-month template (ANTS12-0Months3T), native dimensions | 469 | 133 x 164 | 8088 | 895 | 144 / 82 |

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

Link to G2: at the adult's measured (= centred) position the dense/combined detectability ratio is 1.111x as an unweighted median over all targets (G2's headline), 1.124x without the medial wall and +1.09 dB area-weighted without it (the G3B convention).

## D_child, D_adult and Delta (dense OPM; intrinsic + brain noise; detectability dB)

| child anatomy | comparator | D_child | D_adult | Delta | homology |
|---|---|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | +1.45 dB [+1.23, +1.70] | +0.85 dB [+0.61, +1.04] | +0.59 dB [+0.50, +0.71] | vertex |
| school-age size (scaled adult) | Neuromag grad | +2.96 dB [+2.62, +3.26] | +2.03 dB [+1.83, +2.23] | +0.80 dB [+0.68, +0.97] | vertex |
| school-age size (scaled adult) | Neuromag mag | +1.77 dB [+1.50, +2.06] | +1.23 dB [+0.97, +1.47] | +0.56 dB [+0.47, +0.66] | vertex |
| 2-year size (scaled adult) | Neuromag combined | +2.09 dB [+1.78, +2.39] | +0.85 dB [+0.63, +1.05] | +1.25 dB [+1.11, +1.37] | vertex |
| 2-year size (scaled adult) | Neuromag grad | +3.80 dB [+3.45, +4.11] | +2.03 dB [+1.80, +2.23] | +1.67 dB [+1.50, +1.85] | vertex |
| 2-year size (scaled adult) | Neuromag mag | +2.35 dB [+2.03, +2.66] | +1.23 dB [+0.97, +1.47] | +1.15 dB [+1.01, +1.27] | vertex |
| 2-year template | Neuromag combined | +1.84 dB [+1.52, +2.19] | +0.85 dB [+0.62, +1.05] | +0.88 dB [+0.55, +1.32] | parcel |
| 2-year template | Neuromag grad | +3.33 dB [+2.90, +3.75] | +2.03 dB [+1.78, +2.24] | +1.28 dB [+0.85, +1.75] | parcel |
| 2-year template | Neuromag mag | +2.23 dB [+1.91, +2.57] | +1.23 dB [+1.00, +1.48] | +0.93 dB [+0.63, +1.28] | parcel |
| 18-month template | Neuromag combined | +1.89 dB [+1.61, +2.21] | +0.85 dB [+0.62, +1.05] | +0.99 dB [+0.58, +1.33] | parcel |
| 18-month template | Neuromag grad | +3.38 dB [+2.93, +3.84] | +2.03 dB [+1.81, +2.24] | +1.42 dB [+0.87, +1.86] | parcel |
| 18-month template | Neuromag mag | +2.24 dB [+1.93, +2.60] | +1.23 dB [+0.97, +1.46] | +0.73 dB [+0.67, +1.16] | parcel |
| 12-month template | Neuromag combined | +2.04 dB [+1.64, +2.43] | +0.85 dB [+0.61, +1.06] | +1.23 dB [+0.73, +1.61] | parcel |
| 12-month template | Neuromag grad | +3.65 dB [+3.21, +4.23] | +2.03 dB [+1.79, +2.24] | +1.31 dB [+1.07, +2.28] | parcel |
| 12-month template | Neuromag mag | +2.45 dB [+1.99, +2.84] | +1.23 dB [+0.98, +1.45] | +1.21 dB [+0.77, +1.55] | parcel |

Other metrics, dense OPM vs Neuromag combined, intrinsic + brain (peak-channel SNR: the best single channel, where Neuromag is ahead in the adult; mean-power SNR; both in dB):

| child anatomy | metric | D_child | D_adult | Delta |
|---|---|---|---|---|
| school-age size (scaled adult) | peak-channel SNR | -0.34 dB | -1.30 dB | +0.84 dB |
| school-age size (scaled adult) | mean-power SNR | +0.71 dB | +0.57 dB | +0.14 dB |
| 2-year size (scaled adult) | peak-channel SNR | +0.57 dB | -1.30 dB | +1.83 dB |
| 2-year size (scaled adult) | mean-power SNR | +1.08 dB | +0.57 dB | +0.53 dB |
| 2-year template | peak-channel SNR | +0.32 dB | -1.30 dB | +1.72 dB |
| 2-year template | mean-power SNR | +1.14 dB | +0.57 dB | +0.55 dB |
| 18-month template | peak-channel SNR | +0.45 dB | -1.30 dB | +1.55 dB |
| 18-month template | mean-power SNR | +1.33 dB | +0.57 dB | +0.67 dB |
| 12-month template | peak-channel SNR | +0.79 dB | -1.30 dB | +1.69 dB |
| 12-month template | mean-power SNR | +1.50 dB | +0.57 dB | +0.86 dB |

OPM standoff per anatomy (median sensing-centre height above the MRI scalp, dense / matched array): adult 7.78 / 7.76 mm; school-age size (scaled adult) 7.00 / 7.00 mm; 2-year size (scaled adult) 7.00 / 7.00 mm; 2-year template 7.00 / 7.00 mm; 18-month template 7.00 / 7.01 mm; 12-month template 7.00 / 7.00 mm. The clearance rule (A-OPM-CLEAR) moves most of the adult's sites outward and almost none of the children's (the adult's BEM head surface lies outside its MRI scalp, the templates' inside), so D_adult carries a larger standoff than D_child; the effect on Delta is bounded in `results/g3b/G3B_standoff_report.md` (`scripts/study_g3b_standoff.py`).

Projected condition (room-field subspace removed):

| child anatomy | comparator | D_child | D_adult | Delta |
|---|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | +1.32 dB [+1.04, +1.62] | +0.46 dB [+0.04, +0.76] | +1.10 dB [+0.96, +1.30] |
| school-age size (scaled adult) | Neuromag grad | +2.46 dB [+2.11, +2.80] | +1.28 dB [+0.88, +1.60] | +1.33 dB [+1.15, +1.53] |
| school-age size (scaled adult) | Neuromag mag | +1.83 dB [+1.54, +2.10] | +0.94 dB [+0.54, +1.25] | +1.03 dB [+0.91, +1.18] |
| 2-year size (scaled adult) | Neuromag combined | +2.01 dB [+1.64, +2.35] | +0.46 dB [+0.07, +0.74] | +1.81 dB [+1.69, +1.96] |
| 2-year size (scaled adult) | Neuromag grad | +3.37 dB [+2.90, +3.74] | +1.28 dB [+0.87, +1.60] | +2.22 dB [+2.05, +2.41] |
| 2-year size (scaled adult) | Neuromag mag | +2.46 dB [+2.08, +2.78] | +0.94 dB [+0.55, +1.24] | +1.67 dB [+1.57, +1.82] |
| 2-year template | Neuromag combined | +1.48 dB [+1.04, +1.87] | +0.46 dB [+0.06, +0.75] | +1.00 dB [+0.66, +1.37] |
| 2-year template | Neuromag grad | +2.68 dB [+2.18, +3.22] | +1.28 dB [+0.88, +1.62] | +1.36 dB [+0.93, +1.70] |
| 2-year template | Neuromag mag | +2.02 dB [+1.57, +2.44] | +0.94 dB [+0.54, +1.26] | +1.00 dB [+0.74, +1.37] |
| 18-month template | Neuromag combined | +1.57 dB [+1.20, +1.96] | +0.46 dB [+0.04, +0.75] | +1.16 dB [+0.74, +1.52] |
| 18-month template | Neuromag grad | +2.79 dB [+2.29, +3.28] | +1.28 dB [+0.84, +1.63] | +1.54 dB [+1.02, +1.83] |
| 18-month template | Neuromag mag | +2.06 dB [+1.61, +2.47] | +0.94 dB [+0.56, +1.26] | +1.06 dB [+0.65, +1.32] |
| 12-month template | Neuromag combined | +1.77 dB [+1.29, +2.26] | +0.46 dB [+0.09, +0.74] | +1.46 dB [+0.95, +1.74] |
| 12-month template | Neuromag grad | +3.13 dB [+2.48, +3.64] | +1.28 dB [+0.88, +1.62] | +1.69 dB [+1.17, +2.23] |
| 12-month template | Neuromag mag | +2.32 dB [+1.73, +2.79] | +0.94 dB [+0.55, +1.27] | +1.31 dB [+0.89, +1.66] |

Matched-site OPM (coverage control), intrinsic + brain:

| child anatomy | comparator | D_child | D_adult | Delta |
|---|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | +0.37 dB [+0.21, +0.57] | -0.13 dB [-0.22, -0.03] | +0.53 dB [+0.41, +0.71] |
| school-age size (scaled adult) | Neuromag grad | +1.63 dB [+1.44, +1.89] | +0.74 dB [+0.59, +0.86] | +0.74 dB [+0.58, +0.97] |
| school-age size (scaled adult) | Neuromag mag | +0.67 dB [+0.47, +0.80] | +0.14 dB [-0.01, +0.29] | +0.47 dB [+0.36, +0.65] |
| 2-year size (scaled adult) | Neuromag combined | +1.11 dB [+0.83, +1.35] | -0.13 dB [-0.24, -0.02] | +1.30 dB [+1.15, +1.45] |
| 2-year size (scaled adult) | Neuromag grad | +2.59 dB [+2.25, +2.87] | +0.74 dB [+0.60, +0.87] | +1.71 dB [+1.52, +1.90] |
| 2-year size (scaled adult) | Neuromag mag | +1.38 dB [+1.11, +1.63] | +0.14 dB [-0.01, +0.28] | +1.20 dB [+1.06, +1.35] |
| 2-year template | Neuromag combined | +0.54 dB [+0.29, +0.82] | -0.13 dB [-0.22, -0.04] | +0.93 dB [+0.61, +1.21] |
| 2-year template | Neuromag grad | +1.94 dB [+1.57, +2.38] | +0.74 dB [+0.60, +0.87] | +1.34 dB [+0.88, +1.73] |
| 2-year template | Neuromag mag | +0.87 dB [+0.58, +1.15] | +0.14 dB [-0.02, +0.30] | +0.86 dB [+0.56, +1.08] |
| 18-month template | Neuromag combined | +0.56 dB [+0.33, +0.82] | -0.13 dB [-0.26, -0.04] | +0.91 dB [+0.48, +1.19] |
| 18-month template | Neuromag grad | +1.98 dB [+1.62, +2.33] | +0.74 dB [+0.60, +0.86] | +1.16 dB [+0.92, +1.70] |
| 18-month template | Neuromag mag | +0.87 dB [+0.61, +1.08] | +0.14 dB [-0.01, +0.30] | +0.83 dB [+0.41, +1.04] |
| 12-month template | Neuromag combined | +0.72 dB [+0.45, +1.01] | -0.13 dB [-0.25, -0.03] | +1.14 dB [+0.63, +1.54] |
| 12-month template | Neuromag grad | +2.30 dB [+1.85, +2.84] | +0.74 dB [+0.60, +0.89] | +1.44 dB [+0.96, +2.40] |
| 12-month template | Neuromag mag | +1.05 dB [+0.74, +1.40] | +0.14 dB [-0.01, +0.28] | +0.93 dB [+0.51, +1.48] |

## What drives Delta: placement and helmet fit (dense OPM vs Neuromag combined, intrinsic + brain)

Each child placement is compared with the adult at the same rule (the adult's counterfactual helmet has factor 1).

| child anatomy | top (primary) | centred | x-centred | top-18mm | back | counterfactual | counterfactual, x-centred |
|---|---|---|---|---|---|---|---|
| school-age size (scaled adult) | +0.59 dB [+0.50, +0.71] | +1.29 dB [+1.16, +1.38] | +0.60 dB [+0.50, +0.70] | +0.58 dB [+0.48, +0.70] | +1.03 dB [+0.87, +1.20] | -0.08 dB [-0.11, -0.04] | -0.08 dB [-0.11, -0.04] |
| 2-year size (scaled adult) | +1.25 dB [+1.11, +1.37] | +1.85 dB [+1.69, +1.98] | +1.24 dB [+1.08, +1.36] | +1.28 dB [+1.13, +1.39] | +1.42 dB [+1.12, +1.65] | -0.20 dB [-0.26, -0.16] | -0.21 dB [-0.26, -0.16] |
| 2-year template | +0.88 dB [+0.55, +1.32] | +2.06 dB [+1.60, +2.35] | +0.81 dB [+0.52, +0.97] | +0.89 dB [+0.51, +1.17] | +1.54 dB [+1.46, +1.79] | +0.43 dB [+0.19, +0.82] | -0.20 dB [-0.35, +0.17] |
| 18-month template | +0.99 dB [+0.58, +1.33] | +2.16 dB [+1.87, +2.29] | +0.76 dB [+0.55, +1.16] | +1.01 dB [+0.60, +1.33] | +1.80 dB [+1.27, +2.01] | +0.30 dB [+0.10, +0.73] | -0.12 dB [-0.30, +0.11] |
| 12-month template | +1.23 dB [+0.73, +1.61] | +2.66 dB [+2.18, +3.01] | +1.21 dB [+0.68, +1.51] | +1.12 dB [+0.59, +1.58] | +2.00 dB [+1.38, +2.37] | +0.24 dB [-0.05, +0.37] | +0.01 dB [-0.26, +0.25] |

Counterfactual Delta by comparator (helmet scaled with the head; the dependence on the comparator points to the SQUID side of the change):

| child anatomy | comparator | counterfactual | counterfactual, x-centred |
|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | -0.08 dB [-0.11, -0.04] | -0.08 dB [-0.11, -0.04] |
| school-age size (scaled adult) | Neuromag grad | -0.22 dB [-0.26, -0.17] | -0.22 dB [-0.27, -0.16] |
| school-age size (scaled adult) | Neuromag mag | -0.03 dB [-0.07, +0.01] | -0.03 dB [-0.06, +0.00] |
| 2-year size (scaled adult) | Neuromag combined | -0.20 dB [-0.26, -0.16] | -0.21 dB [-0.26, -0.16] |
| 2-year size (scaled adult) | Neuromag grad | -0.42 dB [-0.49, -0.35] | -0.42 dB [-0.49, -0.35] |
| 2-year size (scaled adult) | Neuromag mag | -0.12 dB [-0.17, -0.08] | -0.13 dB [-0.17, -0.09] |
| 2-year template | Neuromag combined | +0.43 dB [+0.19, +0.82] | -0.20 dB [-0.35, +0.17] |
| 2-year template | Neuromag grad | +0.39 dB [+0.15, +0.62] | -0.26 dB [-0.52, -0.06] |
| 2-year template | Neuromag mag | +0.56 dB [+0.40, +0.83] | +0.00 dB [-0.25, +0.32] |
| 18-month template | Neuromag combined | +0.30 dB [+0.10, +0.73] | -0.12 dB [-0.30, +0.11] |
| 18-month template | Neuromag grad | +0.46 dB [+0.22, +0.60] | -0.14 dB [-0.40, +0.11] |
| 18-month template | Neuromag mag | +0.44 dB [+0.27, +0.69] | -0.07 dB [-0.16, +0.26] |
| 12-month template | Neuromag combined | +0.24 dB [-0.05, +0.37] | +0.01 dB [-0.26, +0.25] |
| 12-month template | Neuromag grad | +0.10 dB [-0.25, +0.24] | -0.13 dB [-0.55, +0.07] |
| 12-month template | Neuromag mag | +0.31 dB [+0.11, +0.49] | +0.12 dB [+0.01, +0.31] |

## Absolute detectability (median 20 log10 d of a 10-nAm dipole, intrinsic + brain, primary placement)

| anatomy | OPM dense | OPM matched | Neuromag combined | Neuromag grad | Neuromag mag |
|---|---|---|---|---|---|
| adult | -1.71 | -2.92 | -2.66 | -3.63 | -3.00 |
| school-age size (scaled adult) | -0.68 | -1.77 | -2.15 | -3.34 | -2.42 |
| 2-year size (scaled adult) | -0.23 | -1.22 | -2.27 | -3.66 | -2.52 |
| 2-year template | +0.92 | -0.41 | -1.19 | -2.70 | -1.52 |
| 18-month template | +1.08 | -0.29 | -1.02 | -2.66 | -1.34 |
| 12-month template | +1.86 | +0.45 | -0.62 | -2.23 | -0.96 |

Vertex-wise change from the adult (scaled controls; same vertex): both systems gain, the OPM more.

| child anatomy | OPM dense | Neuromag combined | Neuromag grad | Neuromag mag |
|---|---|---|---|---|
| school-age size (scaled adult) | +1.25 | +0.65 | +0.45 | +0.69 |
| 2-year size (scaled adult) | +1.83 | +0.54 | +0.10 | +0.64 |

## Channel count: the adult's dense array subsampled to each child's site count

| child anatomy | sites | D_child | D_adult, subsampled | Delta at equal channel count |
|---|---|---|---|---|
| school-age size (scaled adult) | 172 | +1.45 dB [+1.27, +1.68] | +0.63 dB [+0.42, +0.83] | +0.77 dB [+0.63, +0.89] |
| 2-year size (scaled adult) | 155 | +2.09 dB [+1.77, +2.43] | +0.46 dB [+0.28, +0.61] | +1.56 dB [+1.34, +1.69] |
| 2-year template | 151 | +1.84 dB [+1.56, +2.19] | +0.41 dB [+0.27, +0.62] | +1.51 dB [+0.96, +1.78] |
| 18-month template | 157 | +1.89 dB [+1.65, +2.20] | +0.49 dB [+0.33, +0.65] | +1.40 dB [+0.90, +1.76] |
| 12-month template | 144 | +2.04 dB [+1.67, +2.42] | +0.34 dB [+0.21, +0.50] | +1.61 dB [+1.35, +2.44] |

## Delta by depth stratum (dense OPM vs Neuromag combined, intrinsic + brain)

| child anatomy | depth [mm] | n child / adult | D_child | D_adult | Delta [95 % CI] |
|---|---|---|---|---|---|
| school-age size (scaled adult) | 0-10 | 0 / 0 | sparse | | |
| school-age size (scaled adult) | 10-15 | 450 / 222 | +3.87 | +3.40 | +0.47 [-0.05, +0.93] |
| school-age size (scaled adult) | 15-20 | 1604 / 1221 | +2.79 | +2.22 | +0.57 [+0.38, +0.80] |
| school-age size (scaled adult) | 20-25 | 1604 / 1572 | +1.73 | +1.39 | +0.34 [+0.20, +0.53] |
| school-age size (scaled adult) | 25-30 | 1144 / 1289 | +1.02 | +0.77 | +0.25 [+0.12, +0.38] |
| school-age size (scaled adult) | 30-40 | 1399 / 1539 | +0.57 | +0.30 | +0.27 [+0.17, +0.38] |
| school-age size (scaled adult) | 40-50 | 633 / 913 | +0.41 | +0.06 | +0.35 [+0.24, +0.47] |
| school-age size (scaled adult) | 50-60 | 115 / 321 | +0.75 | -0.11 | +0.85 [-0.13, +1.01] |
| school-age size (scaled adult) | 60-90 | 0 / 26 | sparse | | |
| 2-year size (scaled adult) | 0-10 | 2 / 0 | sparse | | |
| 2-year size (scaled adult) | 10-15 | 610 / 222 | +4.88 | +3.40 | +1.48 [+1.00, +1.93] |
| 2-year size (scaled adult) | 15-20 | 1768 / 1221 | +3.44 | +2.22 | +1.22 [+1.01, +1.46] |
| 2-year size (scaled adult) | 20-25 | 1610 / 1572 | +2.13 | +1.39 | +0.74 [+0.59, +0.93] |
| 2-year size (scaled adult) | 25-30 | 1055 / 1289 | +1.24 | +0.77 | +0.47 [+0.31, +0.64] |
| 2-year size (scaled adult) | 30-40 | 1301 / 1539 | +0.74 | +0.30 | +0.44 [+0.33, +0.59] |
| 2-year size (scaled adult) | 40-50 | 463 / 913 | +0.59 | +0.06 | +0.53 [+0.41, +0.67] |
| 2-year size (scaled adult) | 50-60 | 46 / 321 | +1.05 | -0.11 | +1.16 [-0.05, +1.28] |
| 2-year size (scaled adult) | 60-90 | 0 / 26 | sparse | | |
| 2-year template | 0-10 | 26 / 0 | sparse | | |
| 2-year template | 10-15 | 1488 / 222 | +4.08 | +3.40 | +0.67 [+0.13, +1.26] |
| 2-year template | 15-20 | 1666 / 1221 | +2.63 | +2.22 | +0.42 [+0.07, +0.79] |
| 2-year template | 20-25 | 1410 / 1572 | +1.65 | +1.39 | +0.27 [+0.05, +0.49] |
| 2-year template | 25-30 | 980 / 1289 | +1.04 | +0.77 | +0.26 [+0.10, +0.45] |
| 2-year template | 30-40 | 1175 / 1539 | +0.73 | +0.30 | +0.43 [+0.28, +0.70] |
| 2-year template | 40-50 | 502 / 913 | +0.54 | +0.06 | +0.49 [+0.38, +0.84] |
| 2-year template | 50-60 | 229 / 321 | +1.19 | -0.11 | +1.30 [+0.51, +1.60] |
| 2-year template | 60-90 | 83 / 26 | +2.19 | -0.35 | +2.55 [+2.05, +2.74] |
| 18-month template | 0-10 | 32 / 0 | sparse | | |
| 18-month template | 10-15 | 744 / 222 | +4.21 | +3.40 | +0.81 [+0.22, +1.47] |
| 18-month template | 15-20 | 1788 / 1221 | +2.94 | +2.22 | +0.72 [+0.34, +1.30] |
| 18-month template | 20-25 | 1539 / 1572 | +1.92 | +1.39 | +0.54 [+0.21, +0.95] |
| 18-month template | 25-30 | 1161 / 1289 | +1.31 | +0.77 | +0.54 [+0.31, +0.76] |
| 18-month template | 30-40 | 1290 / 1539 | +0.85 | +0.30 | +0.55 [+0.38, +0.78] |
| 18-month template | 40-50 | 523 / 913 | +0.54 | +0.06 | +0.48 [+0.38, +0.79] |
| 18-month template | 50-60 | 221 / 321 | +1.24 | -0.11 | +1.35 [+0.60, +1.62] |
| 18-month template | 60-90 | 65 / 26 | +2.10 | -0.35 | +2.45 [+1.96, +2.57] |
| 12-month template | 0-10 | 78 / 0 | sparse | | |
| 12-month template | 10-15 | 1744 / 222 | +4.55 | +3.40 | +1.15 [+0.41, +1.87] |
| 12-month template | 15-20 | 1768 / 1221 | +2.90 | +2.22 | +0.69 [+0.14, +1.09] |
| 12-month template | 20-25 | 1452 / 1572 | +1.78 | +1.39 | +0.40 [+0.06, +0.71] |
| 12-month template | 25-30 | 928 / 1289 | +1.09 | +0.77 | +0.31 [+0.15, +0.52] |
| 12-month template | 30-40 | 1029 / 1539 | +0.71 | +0.30 | +0.41 [+0.26, +0.60] |
| 12-month template | 40-50 | 390 / 913 | +0.58 | +0.06 | +0.52 [+0.40, +0.72] |
| 12-month template | 50-60 | 184 / 321 | +1.12 | -0.11 | +1.23 [+0.70, +1.38] |
| 12-month template | 60-90 | 19 / 26 | +1.42 | -0.35 | +1.77 [+1.28, +1.90] |

Scaled controls, vertex-wise (homologous) Delta by the adult's depth:

| child anatomy | adult depth [mm] | n | Delta [95 % CI] |
|---|---|---|---|
| school-age size (scaled adult) | 0-10 | 0 | sparse |
| school-age size (scaled adult) | 10-15 | 175 | +1.00 [+0.65, +1.51] |
| school-age size (scaled adult) | 15-20 | 1134 | +0.99 [+0.85, +1.19] |
| school-age size (scaled adult) | 20-25 | 1564 | +0.82 [+0.70, +0.96] |
| school-age size (scaled adult) | 25-30 | 1284 | +0.60 [+0.52, +0.72] |
| school-age size (scaled adult) | 30-40 | 1533 | +0.38 [+0.31, +0.46] |
| school-age size (scaled adult) | 40-50 | 912 | +0.33 [+0.22, +0.49] |
| school-age size (scaled adult) | 50-60 | 321 | +0.53 [+0.29, +0.78] |
| school-age size (scaled adult) | 60-90 | 26 | +1.28 [-0.49, +1.51] |
| 2-year size (scaled adult) | 0-10 | 0 | sparse |
| 2-year size (scaled adult) | 10-15 | 152 | +2.45 [+1.94, +2.67] |
| 2-year size (scaled adult) | 15-20 | 1082 | +2.04 [+1.93, +2.19] |
| 2-year size (scaled adult) | 20-25 | 1551 | +1.68 [+1.54, +1.81] |
| 2-year size (scaled adult) | 25-30 | 1283 | +1.23 [+1.11, +1.34] |
| 2-year size (scaled adult) | 30-40 | 1532 | +0.75 [+0.68, +0.86] |
| 2-year size (scaled adult) | 40-50 | 909 | +0.54 [+0.42, +0.68] |
| 2-year size (scaled adult) | 50-60 | 320 | +0.63 [+0.45, +0.93] |
| 2-year size (scaled adult) | 60-90 | 26 | +1.44 [-0.20, +1.73] |

Templates: the pooled difference of the medians (D_child - D_adult over all targets) and the same with the template's targets reweighted to the adult's area share per depth stratum; then radial (0-30 deg) and tangential (60-90 deg) sources at matched depth (difference of the medians; '-': fewer than 10 targets):

| template | pooled | depth-reweighted | area at 10-20 mm (adult) | median depth [mm] (adult) | radial 0-15 / 15-25 / 25-40 / 40-90 mm | tangential 0-15 / 15-25 / 25-40 / 40-90 mm |
|---|---|---|---|---|---|---|
| 2-year template | +1.00 | +0.48 | 42% (21%) | 21.8 (26.2) | +1.00 / +0.94 / +1.03 / +0.80 | +0.50 / +0.35 / +0.19 / +0.78 |
| 18-month template | +1.05 | +0.68 | 36% (21%) | 23.2 (26.2) | +1.43 / +1.61 / +1.28 / +0.83 | +0.66 / +0.55 / +0.29 / +0.70 |
| 12-month template | +1.20 | +0.44 | 46% (21%) | 20.6 (26.2) | +1.73 / +0.99 / +1.04 / +0.81 | +0.80 / +0.47 / +0.18 / +0.75 |

## Delta by orientation stratum (0 deg = radial to the inner skull; dense OPM vs Neuromag combined)

| child anatomy | orientation [deg] | n child / adult | D_child | D_adult | Delta [95 % CI] |
|---|---|---|---|---|---|
| school-age size (scaled adult) | 0-30 | 676 / 707 | +1.19 | +0.30 | +0.88 [+0.66, +1.18] |
| school-age size (scaled adult) | 30-60 | 2297 / 2346 | +1.23 | +0.61 | +0.62 [+0.30, +0.89] |
| school-age size (scaled adult) | 60-90.1 | 3976 / 4050 | +1.64 | +1.05 | +0.59 [+0.26, +0.94] |
| 2-year size (scaled adult) | 0-30 | 665 / 707 | +1.69 | +0.30 | +1.39 [+1.05, +1.64] |
| 2-year size (scaled adult) | 30-60 | 2269 / 2346 | +1.80 | +0.61 | +1.19 [+0.80, +1.55] |
| 2-year size (scaled adult) | 60-90.1 | 3921 / 4050 | +2.38 | +1.05 | +1.33 [+0.89, +1.71] |
| 2-year template | 0-30 | 1219 / 707 | +2.78 | +0.30 | +2.48 [+1.82, +3.06] |
| 2-year template | 30-60 | 2293 / 2346 | +1.76 | +0.61 | +1.15 [+0.76, +1.51] |
| 2-year template | 60-90.1 | 4047 / 4050 | +1.72 | +1.05 | +0.67 [+0.32, +1.02] |
| 18-month template | 0-30 | 1026 / 707 | +2.57 | +0.30 | +2.26 [+1.68, +2.73] |
| 18-month template | 30-60 | 2398 / 2346 | +1.91 | +0.61 | +1.30 [+0.86, +1.69] |
| 18-month template | 60-90.1 | 3939 / 4050 | +1.77 | +1.05 | +0.72 [+0.41, +1.11] |
| 12-month template | 0-30 | 1294 / 707 | +3.16 | +0.30 | +2.86 [+1.98, +3.86] |
| 12-month template | 30-60 | 2461 / 2346 | +1.90 | +0.61 | +1.29 [+0.79, +1.76] |
| 12-month template | 60-90.1 | 3837 / 4050 | +1.87 | +1.05 | +0.82 [+0.42, +1.24] |

## Placement, counterfactual helmet and sensitivity (median D, dense OPM vs Neuromag combined, intrinsic + brain)

| anatomy | centred | top | back | x+5mm | x-5mm | y+5mm | y-5mm | pitch+10deg | pitch-10deg | roll+5deg | roll-5deg | x-centred | top-18mm | counterfactual | counterfactual_x-centred |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| adult | +1.09 | +0.85 | +0.98 | +0.96 | +1.12 | +1.09 | +0.83 | +0.80 | +1.25 | +1.10 | +1.00 | +0.87 | +0.72 | +1.09 | +1.11 |
| school-age size (scaled adult) | +2.43 | +1.45 | +1.94 | +1.44 | +1.49 | +1.52 | +1.41 | +1.45 | +2.48 | +1.42 | +1.65 | +1.46 | +1.32 | +0.90 | +0.91 |
| 2-year size (scaled adult) | +3.00 | +2.09 | +2.29 | +1.96 | +2.11 | +2.17 | +2.01 | +1.64 | +2.98 | +1.88 | +2.45 | +2.09 | +1.98 | +0.68 | +0.69 |
| 2-year template | +3.08 | +1.84 | +2.51 | +2.07 | +1.75 | +1.84 | +1.85 | +1.78 | +2.00 | +1.79 | +2.04 | +1.76 | +1.70 | +1.55 | +1.00 |
| 18-month template | +3.13 | +1.89 | +2.48 | +2.05 | +1.80 | +1.86 | +1.82 | +1.80 | +2.03 | +1.80 | +2.04 | +1.81 | +1.75 | +1.47 | +1.02 |
| 12-month template | +3.65 | +2.04 | +2.79 | +2.13 | +2.03 | +2.07 | +2.06 | +1.98 | +2.22 | +2.01 | +2.11 | +1.98 | +1.85 | +1.16 | +0.99 |

| anatomy | opm_asd_7fT | opm_asd_10fT | opm_asd_15fT | opm_asd_20fT | opm_asd_30fT | background_x0.5 | background_x2 | bem1 |
|---|---|---|---|---|---|---|---|---|
| adult | +1.80 | +1.37 | +0.85 | +0.44 | -0.15 | +0.82 | +0.85 | +0.90 |
| school-age size (scaled adult) | +2.48 | +2.05 | +1.45 | +0.99 | +0.31 | +1.44 | +1.44 | +1.43 |
| 2-year size (scaled adult) | +3.11 | +2.66 | +2.09 | +1.64 | +0.94 | +2.10 | +2.03 | +2.08 |
| 2-year template | +2.81 | +2.40 | +1.84 | +1.38 | +0.61 | +1.86 | +1.79 | +1.80 |
| 18-month template | +2.98 | +2.53 | +1.89 | +1.39 | +0.61 | +1.86 | +1.90 | +1.81 |
| 12-month template | +3.06 | +2.64 | +2.04 | +1.54 | +0.77 | +2.05 | +1.98 | +2.03 |

Difference of these medians, child minus adult, with the same variant applied to both (a sensitivity of the medians, not the paired Delta estimator; the background variants scale the adult too):

| child anatomy | opm_asd_7fT | opm_asd_10fT | opm_asd_15fT | opm_asd_20fT | opm_asd_30fT | background_x0.5 | background_x2 | bem1 |
|---|---|---|---|---|---|---|---|---|
| school-age size (scaled adult) | +0.67 | +0.68 | +0.61 | +0.55 | +0.46 | +0.62 | +0.59 | +0.53 |
| 2-year size (scaled adult) | +1.31 | +1.29 | +1.24 | +1.20 | +1.09 | +1.28 | +1.18 | +1.17 |
| 2-year template | +1.01 | +1.03 | +1.00 | +0.94 | +0.76 | +1.04 | +0.94 | +0.90 |
| 18-month template | +1.18 | +1.16 | +1.05 | +0.95 | +0.76 | +1.04 | +1.05 | +0.91 |
| 12-month template | +1.26 | +1.27 | +1.20 | +1.10 | +0.92 | +1.23 | +1.13 | +1.13 |

## Regions that gain or lose with the placement (median D by lobe, dense OPM vs Neuromag combined, intrinsic + brain)

| anatomy | placement | frontal | parietal | temporal | occipital | cingulate | insula |
|---|---|---|---|---|---|---|---|
| adult | centred | +1.71 | +1.05 | +0.92 | +0.65 | +0.13 | +0.31 |
| adult | top | +1.33 | +0.75 | +0.80 | +0.54 | +0.02 | +0.13 |
| adult | x-centred | +1.37 | +0.77 | +0.80 | +0.54 | +0.02 | +0.10 |
| adult | back | +1.85 | +0.85 | +0.78 | +0.30 | +0.11 | +0.22 |
| adult | counterfactual | +1.71 | +1.05 | +0.92 | +0.65 | +0.13 | +0.31 |
| school-age size (scaled adult) | centred | +3.06 | +2.53 | +2.38 | +1.98 | +0.87 | +1.24 |
| school-age size (scaled adult) | top | +1.71 | +1.33 | +1.77 | +1.49 | +0.35 | +0.76 |
| school-age size (scaled adult) | x-centred | +1.72 | +1.33 | +1.78 | +1.49 | +0.35 | +0.76 |
| school-age size (scaled adult) | back | +3.36 | +1.86 | +1.86 | +0.72 | +0.74 | +1.06 |
| school-age size (scaled adult) | counterfactual | +1.33 | +0.90 | +0.81 | +0.59 | +0.28 | +0.44 |
| 2-year size (scaled adult) | centred | +3.73 | +3.05 | +3.03 | +2.35 | +1.07 | +1.60 |
| 2-year size (scaled adult) | top | +2.55 | +2.02 | +2.41 | +1.90 | +0.56 | +1.12 |
| 2-year size (scaled adult) | x-centred | +2.56 | +2.02 | +2.40 | +1.89 | +0.56 | +1.13 |
| 2-year size (scaled adult) | back | +3.99 | +2.11 | +2.21 | +0.74 | +0.86 | +1.31 |
| 2-year size (scaled adult) | counterfactual | +1.14 | +0.66 | +0.59 | +0.36 | +0.17 | +0.34 |
| 2-year template | centred | +3.71 | +3.02 | +3.27 | +2.74 | +1.10 | +1.96 |
| 2-year template | top | +2.02 | +1.39 | +2.63 | +2.08 | +0.61 | +1.34 |
| 2-year template | x-centred | +1.95 | +1.30 | +2.63 | +2.02 | +0.57 | +1.28 |
| 2-year template | back | +3.89 | +2.20 | +2.77 | +1.25 | +1.07 | +1.74 |
| 2-year template | counterfactual | +2.16 | +1.25 | +1.72 | +1.21 | +0.66 | +1.15 |
| 18-month template | centred | +3.86 | +2.80 | +3.51 | +2.58 | +1.07 | +1.88 |
| 18-month template | top | +2.25 | +1.34 | +2.88 | +1.91 | +0.56 | +1.37 |
| 18-month template | x-centred | +2.15 | +1.23 | +2.84 | +1.85 | +0.53 | +1.27 |
| 18-month template | back | +4.03 | +1.99 | +2.80 | +1.20 | +1.04 | +1.75 |
| 18-month template | counterfactual | +2.14 | +1.03 | +1.76 | +1.00 | +0.57 | +1.07 |
| 12-month template | centred | +4.51 | +3.59 | +3.80 | +3.20 | +1.18 | +1.81 |
| 12-month template | top | +2.46 | +1.44 | +3.02 | +2.24 | +0.59 | +1.24 |
| 12-month template | x-centred | +2.42 | +1.37 | +3.01 | +2.26 | +0.57 | +1.19 |
| 12-month template | back | +4.69 | +2.58 | +2.92 | +1.07 | +1.01 | +1.59 |
| 12-month template | counterfactual | +2.06 | +0.85 | +1.24 | +0.62 | +0.44 | +0.70 |

## Source-to-sensor distance (median, mm) by depth below the scalp

| anatomy | sensors | 0-10 | 10-15 | 15-20 | 20-25 | 25-30 | 30-40 | 40-50 | 50-60 | 60-90 |
|---|---|---|---|---|---|---|---|---|---|---|
| adult | squid:top | - | 45 | 48 | 52 | 57 | 64 | 74 | 82 | 87 |
| adult | opm_dense | - | 23 | 27 | 32 | 36 | 44 | 54 | 64 | 69 |
| school-age size (scaled adult) | squid:top | - | 48 | 52 | 56 | 61 | 69 | 75 | 78 | - |
| school-age size (scaled adult) | opm_dense | - | 22 | 26 | 31 | 36 | 43 | 53 | 61 | - |
| 2-year size (scaled adult) | squid:top | 57 | 54 | 57 | 62 | 67 | 75 | 80 | 84 | - |
| 2-year size (scaled adult) | opm_dense | 18 | 22 | 26 | 31 | 36 | 43 | 53 | 59 | - |
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
| adult | opm_dense/mag | 490.9 | 88.8 |
| adult | opm_matched/mag | 505.0 | 88.8 |
| school-age size (scaled adult) | squid:top/grad | 29.8 | 21.3 |
| school-age size (scaled adult) | squid:top/mag | 149.8 | 20.7 |
| school-age size (scaled adult) | squid:centred/grad | 21.4 | 21.3 |
| school-age size (scaled adult) | squid:centred/mag | 131.0 | 20.7 |
| school-age size (scaled adult) | opm_dense/mag | 515.7 | 88.8 |
| school-age size (scaled adult) | opm_matched/mag | 515.0 | 88.8 |
| 2-year size (scaled adult) | squid:top/grad | 23.4 | 21.3 |
| 2-year size (scaled adult) | squid:top/mag | 124.9 | 20.7 |
| 2-year size (scaled adult) | squid:centred/grad | 16.6 | 21.3 |
| 2-year size (scaled adult) | squid:centred/mag | 108.0 | 20.7 |
| 2-year size (scaled adult) | opm_dense/mag | 528.3 | 88.8 |
| 2-year size (scaled adult) | opm_matched/mag | 541.7 | 88.8 |
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
| 12-month template | opm_dense/mag | 474.1 | 88.8 |
| 12-month template | opm_matched/mag | 470.0 | 88.8 |

## Usefulness (dense OPM vs Neuromag combined, intrinsic + brain): share of usable cortical area

A source counts as usable when its detectability reaches 5 at the reference moment (an operational choice, not a clinical standard).

| anatomy | moment [nAm] | both | OPM only | SQUID only | neither |
|---|---|---|---|---|---|
| adult | 20 | 0.00 | 0.03 | 0.00 | 0.97 |
| adult | 50 | 0.35 | 0.07 | 0.00 | 0.58 |
| adult | 100 | 0.66 | 0.02 | 0.00 | 0.32 |
| adult | 200 | 0.89 | 0.00 | 0.00 | 0.11 |
| school-age size (scaled adult) | 20 | 0.00 | 0.06 | 0.00 | 0.94 |
| school-age size (scaled adult) | 50 | 0.37 | 0.10 | 0.00 | 0.53 |
| school-age size (scaled adult) | 100 | 0.69 | 0.03 | 0.00 | 0.28 |
| school-age size (scaled adult) | 200 | 0.91 | 0.02 | 0.00 | 0.07 |
| 2-year size (scaled adult) | 20 | 0.00 | 0.08 | 0.00 | 0.92 |
| 2-year size (scaled adult) | 50 | 0.35 | 0.14 | 0.00 | 0.51 |
| 2-year size (scaled adult) | 100 | 0.69 | 0.05 | 0.00 | 0.26 |
| 2-year size (scaled adult) | 200 | 0.91 | 0.02 | 0.00 | 0.07 |
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
- Intervals: bootstrap over parcels of one anatomy (or of each anatomy, for between-anatomy strata); they do not include between-subject variability. Three average templates of one database (2-year template, 18-month template, 12-month template) are not a population: template results are conditional simulations.
- Every child array uses the adult's conventions: background moment variance per unit cortical area, room field, intrinsic noise, sensor sizes and the 3-layer BEM conductivities; only geometry changes. Both systems' detectability rises in the smaller heads, the OPM's more (absolute detectability table), by different routes: the on-scalp OPM sees more signal from a cortex that is closer in absolute terms at about the same brain noise, while the SQUIDs' brain noise falls (the cortex is farther from the fixed helmet and, with the background fixed per unit area, smaller) more than their signal. The templates' averaged white surfaces are smoother than an individual cortex (usable area 1,062, 975, 895 cm^2 for the 2-year template, 18-month template, 12-month template vs 1,878 cm^2 for the adult), which lowers their background power and their patch cancellation further; scaling the background variance x0.5 or x2 leaves D_child almost unchanged.
- Placements are chosen from the scalp and helmet geometry only. Under the adult's measured pose a head with other fiducials need not be centred laterally; 'x-centred' shifts each head along device x to equal left/right median gaps before the top contact (shift: adult -1.5 mm, school-age size (scaled adult) -0.5 mm, 2-year size (scaled adult) -0.5 mm, 2-year template -6.5 mm, 18-month template -5.5 mm, 12-month template -2.5 mm; negative = to the left), and 'counterfactual_x-centred' scales the helmet about that laterally centred head. The counterfactual helmet (scaled with the head) is a mechanistic control, not a pediatric SQUID system.
- Targets on the medial wall (FreeSurfer 'unknown': the cut through the corpus callosum and midbrain, not cortex) are left out of every summary; they would otherwise dominate the deepest strata.
