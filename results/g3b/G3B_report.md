# G3B: fixed adult Neuromag helmet versus head-adaptive OPM on smaller heads (NEW)

Code commit: see `g3b_summary.json` (provenance). Configuration: `configs/g3b_pediatric.toml`. Primary placement: `top` (head raised to 20-mm contact). Metric: known-topography detectability in dB; D = OPM - SQUID; Delta = D_child - D_adult. Intervals: parcel bootstrap (one anatomy each; no between-subject variability).

## Anatomies

| anatomy | description | OFC [mm] | breadth x length [mm] | targets | usable cortex [cm2] | OPM dense / matched sites |
|---|---|---|---|---|---|---|
| adult | MNE sample subject (as G2) | 586 | 171 x 211 | 7661 | 1878 | 208 / 98 |
| school-age size (scaled adult) | adult x 85/95 (Jas Table 1 child/adult head radius): 0.8947 | 525 | 153 x 189 | 7507 | 1470 | 174 / 89 |
| 2-year size (scaled adult) | adult x template/adult occipitofrontal circumference: 0.8442 | 495 | 144 x 178 | 7413 | 1290 | 155 / 90 |
| 2-year template | 2-year template (ANTS2-0Years3T), native dimensions | 495 | 140 x 175 | 8134 | 1062 | 151 / 83 |
| 18-month template | 18-month template (ANTS18-0Months3T), native dimensions | 491 | 137 x 172 | 7994 | 975 | 157 / 80 |
| 12-month template | 12-month template (ANTS12-0Months3T), native dimensions | 469 | 133 x 164 | 8088 | 895 | 144 / 82 |
| child A (7.8 y) | child A (7.8 y) (sub-Z213), individual MRI, modelled skull | 520 | 154 x 177 | 7618 | 1852 | 155 / 90 |
| child B (8.3 y) | child B (8.3 y) (sub-Z209), individual MRI, modelled skull | 486 | 141 x 174 | 7428 | 1666 | 153 / 85 |
| child C (8.7 y) | child C (8.7 y) (sub-Z226), individual MRI, modelled skull | 535 | 161 x 185 | 7947 | 1869 | 167 / 91 |

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
| child A (7.8 y) | centred | 0.0 | 30.9 | 40.2 | True |
| child A (7.8 y) | top | 17.0 | 20.1 | 35.2 | True |
| child A (7.8 y) | back | 21.0 | 20.1 | 36.0 | True |
| child A (7.8 y) | x-centred | 17.0 | 20.1 | 35.2 | True |
| child A (7.8 y) | top-18mm | 19.0 | 18.4 | 34.8 | True |
| child A (7.8 y) | counterfactual | 0.0 | 18.1 | 26.1 | True |
| child A (7.8 y) | counterfactual_x-centred | 0.0 | 18.1 | 26.1 | True |
| child B (8.3 y) | centred | 0.0 | 33.6 | 47.3 | True |
| child B (8.3 y) | top | 26.5 | 20.1 | 38.4 | True |
| child B (8.3 y) | back | 24.5 | 20.1 | 42.2 | True |
| child B (8.3 y) | x-centred | 27.5 | 20.3 | 38.7 | True |
| child B (8.3 y) | top-18mm | 28.5 | 18.3 | 38.1 | True |
| child B (8.3 y) | counterfactual | 0.0 | 18.5 | 31.4 | True |
| child B (8.3 y) | counterfactual_x-centred | 0.0 | 18.1 | 28.0 | True |
| child C (8.7 y) | centred | 0.0 | 21.6 | 41.3 | True |
| child C (8.7 y) | top | 21.0 | 20.0 | 32.6 | True |
| child C (8.7 y) | back | 11.0 | 20.0 | 38.0 | True |
| child C (8.7 y) | x-5mm | 0.0 | 16.9 | 42.5 | False |
| child C (8.7 y) | x-centred | 20.5 | 20.2 | 33.1 | True |
| child C (8.7 y) | top-18mm | 23.0 | 18.1 | 31.9 | True |
| child C (8.7 y) | counterfactual | 0.0 | 18.5 | 36.8 | True |
| child C (8.7 y) | counterfactual_x-centred | 0.0 | 18.2 | 33.1 | True |

Link to G2: at the adult's measured (= centred) position the dense/combined detectability ratio is 1.144x as an unweighted median over all targets (G2's headline), 1.148x without the medial wall and +1.27 dB area-weighted without it (the G3B convention).

## D_child, D_adult and Delta (dense OPM; intrinsic + brain noise; detectability dB)

| child anatomy | comparator | D_child | D_adult | Delta | homology |
|---|---|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | +1.50 dB [+1.28, +1.75] | +1.00 dB [+0.82, +1.17] | +0.44 dB [+0.34, +0.56] | vertex |
| school-age size (scaled adult) | Neuromag grad | +3.00 dB [+2.68, +3.30] | +2.23 dB [+2.03, +2.45] | +0.65 dB [+0.51, +0.80] | vertex |
| school-age size (scaled adult) | Neuromag mag | +1.81 dB [+1.55, +2.12] | +1.39 dB [+1.18, +1.61] | +0.41 dB [+0.32, +0.51] | vertex |
| 2-year size (scaled adult) | Neuromag combined | +2.13 dB [+1.83, +2.44] | +1.00 dB [+0.84, +1.18] | +1.07 dB [+0.91, +1.18] | vertex |
| 2-year size (scaled adult) | Neuromag grad | +3.83 dB [+3.48, +4.14] | +2.23 dB [+2.02, +2.44] | +1.49 dB [+1.31, +1.67] | vertex |
| 2-year size (scaled adult) | Neuromag mag | +2.41 dB [+2.07, +2.71] | +1.39 dB [+1.17, +1.60] | +0.96 dB [+0.82, +1.10] | vertex |
| 2-year template | Neuromag combined | +1.84 dB [+1.52, +2.19] | +1.00 dB [+0.83, +1.18] | +0.73 dB [+0.46, +0.98] | parcel |
| 2-year template | Neuromag grad | +3.33 dB [+2.90, +3.75] | +2.23 dB [+2.00, +2.45] | +1.13 dB [+0.77, +1.55] | parcel |
| 2-year template | Neuromag mag | +2.23 dB [+1.91, +2.57] | +1.39 dB [+1.20, +1.61] | +0.82 dB [+0.55, +0.93] | parcel |
| 18-month template | Neuromag combined | +1.89 dB [+1.61, +2.21] | +1.00 dB [+0.83, +1.18] | +0.88 dB [+0.54, +1.05] | parcel |
| 18-month template | Neuromag grad | +3.38 dB [+2.93, +3.84] | +2.23 dB [+2.02, +2.45] | +1.26 dB [+0.77, +1.75] | parcel |
| 18-month template | Neuromag mag | +2.24 dB [+1.93, +2.60] | +1.39 dB [+1.16, +1.61] | +0.72 dB [+0.55, +1.04] | parcel |
| 12-month template | Neuromag combined | +2.04 dB [+1.64, +2.43] | +1.00 dB [+0.83, +1.18] | +0.96 dB [+0.61, +1.51] | parcel |
| 12-month template | Neuromag grad | +3.65 dB [+3.21, +4.23] | +2.23 dB [+2.01, +2.45] | +1.33 dB [+0.82, +2.18] | parcel |
| 12-month template | Neuromag mag | +2.45 dB [+2.00, +2.84] | +1.39 dB [+1.19, +1.59] | +0.98 dB [+0.68, +1.47] | parcel |
| child A (7.8 y) | Neuromag combined | +1.34 dB [+1.09, +1.59] | +1.00 dB [+0.83, +1.19] | +0.30 dB [+0.11, +0.38] | parcel |
| child A (7.8 y) | Neuromag grad | +2.34 dB [+2.02, +2.72] | +2.23 dB [+2.04, +2.45] | +0.18 dB [+0.04, +0.28] | parcel |
| child A (7.8 y) | Neuromag mag | +1.71 dB [+1.44, +1.99] | +1.39 dB [+1.19, +1.59] | +0.33 dB [+0.14, +0.37] | parcel |
| child B (8.3 y) | Neuromag combined | +1.51 dB [+1.23, +1.79] | +1.00 dB [+0.82, +1.18] | +0.41 dB [+0.28, +0.54] | parcel |
| child B (8.3 y) | Neuromag grad | +2.59 dB [+2.18, +2.95] | +2.23 dB [+2.03, +2.45] | +0.34 dB [+0.26, +0.53] | parcel |
| child B (8.3 y) | Neuromag mag | +1.93 dB [+1.62, +2.21] | +1.39 dB [+1.18, +1.59] | +0.50 dB [+0.29, +0.66] | parcel |
| child C (8.7 y) | Neuromag combined | +1.19 dB [+0.96, +1.39] | +1.00 dB [+0.82, +1.18] | +0.17 dB [-0.07, +0.34] | parcel |
| child C (8.7 y) | Neuromag grad | +2.30 dB [+1.99, +2.64] | +2.23 dB [+2.03, +2.45] | +0.11 dB [-0.00, +0.32] | parcel |
| child C (8.7 y) | Neuromag mag | +1.55 dB [+1.26, +1.79] | +1.39 dB [+1.17, +1.61] | +0.25 dB [-0.14, +0.37] | parcel |

Other metrics, dense OPM vs Neuromag combined, intrinsic + brain (peak-channel SNR: the best single channel, where Neuromag is ahead in the adult; mean-power SNR; both in dB):

| child anatomy | metric | D_child | D_adult | Delta |
|---|---|---|---|---|
| school-age size (scaled adult) | peak-channel SNR | -0.29 dB | -1.11 dB | +0.69 dB |
| school-age size (scaled adult) | mean-power SNR | +0.74 dB | +0.34 dB | +0.39 dB |
| 2-year size (scaled adult) | peak-channel SNR | +0.59 dB | -1.11 dB | +1.66 dB |
| 2-year size (scaled adult) | mean-power SNR | +1.09 dB | +0.34 dB | +0.74 dB |
| 2-year template | peak-channel SNR | +0.32 dB | -1.11 dB | +1.52 dB |
| 2-year template | mean-power SNR | +1.14 dB | +0.34 dB | +0.83 dB |
| 18-month template | peak-channel SNR | +0.45 dB | -1.11 dB | +1.32 dB |
| 18-month template | mean-power SNR | +1.33 dB | +0.34 dB | +0.92 dB |
| 12-month template | peak-channel SNR | +0.79 dB | -1.11 dB | +1.60 dB |
| 12-month template | mean-power SNR | +1.50 dB | +0.34 dB | +1.08 dB |
| child A (7.8 y) | peak-channel SNR | -0.64 dB | -1.11 dB | +0.64 dB |
| child A (7.8 y) | mean-power SNR | +0.49 dB | +0.34 dB | +0.16 dB |
| child B (8.3 y) | peak-channel SNR | -0.28 dB | -1.11 dB | +0.98 dB |
| child B (8.3 y) | mean-power SNR | +0.47 dB | +0.34 dB | +0.23 dB |
| child C (8.7 y) | peak-channel SNR | -0.85 dB | -1.11 dB | +0.24 dB |
| child C (8.7 y) | mean-power SNR | +0.41 dB | +0.34 dB | +0.17 dB |

OPM standoff per anatomy (median sensing-centre height above the MRI scalp, dense / matched array, and the sites the clearance rule moved outward): adult 7.00 / 6.99 mm (26 of 208 / 9 of 98 moved; exact cell clearance >= 1.00 mm); school-age size (scaled adult) 6.99 / 6.99 mm (18 of 174 / 9 of 89 moved; exact cell clearance >= 1.00 mm); 2-year size (scaled adult) 6.99 / 6.99 mm (17 of 155 / 11 of 90 moved; exact cell clearance >= 1.01 mm); 2-year template 7.00 / 7.00 mm (3 of 151 / 3 of 83 moved; exact cell clearance >= 1.00 mm); 18-month template 7.00 / 7.01 mm (5 of 157 / 4 of 80 moved; exact cell clearance >= 1.00 mm); 12-month template 7.00 / 7.00 mm (4 of 144 / 5 of 82 moved; exact cell clearance >= 1.01 mm); child A (7.8 y) 7.00 / 7.00 mm (6 of 155 / 4 of 90 moved; exact cell clearance >= 0.99 mm); child B (8.3 y) 7.00 / 7.00 mm (3 of 153 / 3 of 85 moved; exact cell clearance >= 1.07 mm); child C (8.7 y) 7.00 / 7.00 mm (4 of 167 / 3 of 91 moved; exact cell clearance >= 1.02 mm). Every anatomy's BEM head surface has its vertices on its MRI scalp (A-BEM-CONFORM: the templates' are built so, the adult's stored outer skin, about 1 mm outside its scalp, is conformed at loading), so the clearance rule (A-OPM-CLEAR) moves only the sites the anatomy demands and every array has the same nominal standoff.

Projected condition (room-field subspace removed):

| child anatomy | comparator | D_child | D_adult | Delta |
|---|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | +1.38 dB [+1.08, +1.68] | +0.86 dB [+0.63, +1.08] | +0.44 dB [+0.32, +0.57] |
| school-age size (scaled adult) | Neuromag grad | +2.54 dB [+2.15, +2.83] | +1.80 dB [+1.56, +2.05] | +0.60 dB [+0.45, +0.76] |
| school-age size (scaled adult) | Neuromag mag | +1.89 dB [+1.59, +2.16] | +1.46 dB [+1.25, +1.68] | +0.34 dB [+0.23, +0.45] |
| 2-year size (scaled adult) | Neuromag combined | +2.03 dB [+1.66, +2.38] | +0.86 dB [+0.66, +1.07] | +1.09 dB [+0.93, +1.23] |
| 2-year size (scaled adult) | Neuromag grad | +3.39 dB [+2.91, +3.78] | +1.80 dB [+1.55, +2.05] | +1.47 dB [+1.26, +1.65] |
| 2-year size (scaled adult) | Neuromag mag | +2.47 dB [+2.08, +2.81] | +1.46 dB [+1.24, +1.68] | +0.94 dB [+0.76, +1.09] |
| 2-year template | Neuromag combined | +1.48 dB [+1.04, +1.87] | +0.86 dB [+0.64, +1.07] | +0.42 dB [+0.15, +0.73] |
| 2-year template | Neuromag grad | +2.68 dB [+2.17, +3.22] | +1.80 dB [+1.55, +2.04] | +0.78 dB [+0.47, +1.17] |
| 2-year template | Neuromag mag | +2.02 dB [+1.57, +2.44] | +1.46 dB [+1.26, +1.68] | +0.47 dB [+0.08, +0.79] |
| 18-month template | Neuromag combined | +1.57 dB [+1.20, +1.96] | +0.86 dB [+0.64, +1.07] | +0.57 dB [+0.31, +0.94] |
| 18-month template | Neuromag grad | +2.79 dB [+2.29, +3.28] | +1.80 dB [+1.54, +2.07] | +1.07 dB [+0.58, +1.60] |
| 18-month template | Neuromag mag | +2.06 dB [+1.61, +2.47] | +1.46 dB [+1.26, +1.69] | +0.52 dB [+0.29, +0.81] |
| 12-month template | Neuromag combined | +1.77 dB [+1.29, +2.26] | +0.86 dB [+0.66, +1.07] | +0.83 dB [+0.35, +1.38] |
| 12-month template | Neuromag grad | +3.13 dB [+2.48, +3.64] | +1.80 dB [+1.55, +2.05] | +1.10 dB [+0.60, +2.01] |
| 12-month template | Neuromag mag | +2.32 dB [+1.73, +2.79] | +1.46 dB [+1.25, +1.69] | +0.83 dB [+0.47, +1.48] |
| child A (7.8 y) | Neuromag combined | +0.95 dB [+0.57, +1.25] | +0.86 dB [+0.64, +1.08] | +0.03 dB [-0.17, +0.22] |
| child A (7.8 y) | Neuromag grad | +1.76 dB [+1.31, +2.14] | +1.80 dB [+1.55, +2.05] | +0.08 dB [-0.21, +0.23] |
| child A (7.8 y) | Neuromag mag | +1.46 dB [+0.99, +1.80] | +1.46 dB [+1.26, +1.69] | +0.03 dB [-0.18, +0.27] |
| child B (8.3 y) | Neuromag combined | +1.20 dB [+0.89, +1.55] | +0.86 dB [+0.66, +1.07] | +0.26 dB [+0.04, +0.39] |
| child B (8.3 y) | Neuromag grad | +2.12 dB [+1.63, +2.56] | +1.80 dB [+1.54, +2.06] | +0.29 dB [-0.00, +0.46] |
| child B (8.3 y) | Neuromag mag | +1.75 dB [+1.33, +2.10] | +1.46 dB [+1.26, +1.68] | +0.28 dB [+0.00, +0.46] |
| child C (8.7 y) | Neuromag combined | +0.90 dB [+0.59, +1.19] | +0.86 dB [+0.66, +1.07] | -0.19 dB [-0.38, +0.18] |
| child C (8.7 y) | Neuromag grad | +1.77 dB [+1.34, +2.14] | +1.80 dB [+1.55, +2.05] | -0.11 dB [-0.31, +0.18] |
| child C (8.7 y) | Neuromag mag | +1.37 dB [+0.98, +1.67] | +1.46 dB [+1.26, +1.68] | -0.22 dB [-0.58, +0.09] |

Matched-site OPM (coverage control), intrinsic + brain:

| child anatomy | comparator | D_child | D_adult | Delta |
|---|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | +0.40 dB [+0.23, +0.59] | -0.01 dB [-0.08, +0.08] | +0.34 dB [+0.25, +0.49] |
| school-age size (scaled adult) | Neuromag grad | +1.64 dB [+1.45, +1.91] | +0.97 dB [+0.79, +1.13] | +0.54 dB [+0.39, +0.74] |
| school-age size (scaled adult) | Neuromag mag | +0.68 dB [+0.49, +0.82] | +0.27 dB [+0.14, +0.41] | +0.31 dB [+0.22, +0.42] |
| 2-year size (scaled adult) | Neuromag combined | +1.12 dB [+0.85, +1.37] | -0.01 dB [-0.08, +0.10] | +1.10 dB [+0.94, +1.25] |
| 2-year size (scaled adult) | Neuromag grad | +2.60 dB [+2.26, +2.87] | +0.97 dB [+0.81, +1.15] | +1.49 dB [+1.28, +1.66] |
| 2-year size (scaled adult) | Neuromag mag | +1.39 dB [+1.13, +1.64] | +0.27 dB [+0.13, +0.40] | +1.00 dB [+0.87, +1.14] |
| 2-year template | Neuromag combined | +0.54 dB [+0.29, +0.82] | -0.01 dB [-0.08, +0.07] | +0.68 dB [+0.46, +0.84] |
| 2-year template | Neuromag grad | +1.94 dB [+1.57, +2.38] | +0.97 dB [+0.83, +1.13] | +1.00 dB [+0.53, +1.48] |
| 2-year template | Neuromag mag | +0.87 dB [+0.58, +1.15] | +0.27 dB [+0.13, +0.41] | +0.68 dB [+0.41, +0.87] |
| 18-month template | Neuromag combined | +0.56 dB [+0.33, +0.82] | -0.01 dB [-0.09, +0.08] | +0.61 dB [+0.29, +0.99] |
| 18-month template | Neuromag grad | +1.98 dB [+1.62, +2.33] | +0.97 dB [+0.81, +1.14] | +1.02 dB [+0.49, +1.51] |
| 18-month template | Neuromag mag | +0.87 dB [+0.61, +1.08] | +0.27 dB [+0.13, +0.42] | +0.60 dB [+0.32, +0.83] |
| 12-month template | Neuromag combined | +0.72 dB [+0.45, +1.01] | -0.01 dB [-0.09, +0.08] | +0.76 dB [+0.44, +1.16] |
| 12-month template | Neuromag grad | +2.30 dB [+1.85, +2.84] | +0.97 dB [+0.81, +1.16] | +1.18 dB [+0.77, +2.08] |
| 12-month template | Neuromag mag | +1.05 dB [+0.74, +1.40] | +0.27 dB [+0.15, +0.40] | +0.66 dB [+0.45, +1.20] |
| child A (7.8 y) | Neuromag combined | +0.27 dB [+0.11, +0.48] | -0.01 dB [-0.08, +0.06] | +0.27 dB [+0.09, +0.36] |
| child A (7.8 y) | Neuromag grad | +1.05 dB [+0.80, +1.40] | +0.97 dB [+0.82, +1.16] | +0.08 dB [-0.08, +0.35] |
| child A (7.8 y) | Neuromag mag | +0.56 dB [+0.37, +0.81] | +0.27 dB [+0.15, +0.40] | +0.25 dB [+0.06, +0.35] |
| child B (8.3 y) | Neuromag combined | +0.21 dB [+0.05, +0.46] | -0.01 dB [-0.09, +0.08] | +0.08 dB [-0.05, +0.31] |
| child B (8.3 y) | Neuromag grad | +1.06 dB [+0.78, +1.40] | +0.97 dB [+0.84, +1.16] | +0.08 dB [-0.26, +0.23] |
| child B (8.3 y) | Neuromag mag | +0.52 dB [+0.32, +0.77] | +0.27 dB [+0.15, +0.40] | +0.11 dB [-0.08, +0.23] |
| child C (8.7 y) | Neuromag combined | +0.17 dB [-0.02, +0.32] | -0.01 dB [-0.08, +0.08] | +0.17 dB [+0.08, +0.29] |
| child C (8.7 y) | Neuromag grad | +1.11 dB [+0.80, +1.37] | +0.97 dB [+0.81, +1.14] | +0.20 dB [-0.07, +0.42] |
| child C (8.7 y) | Neuromag mag | +0.47 dB [+0.26, +0.62] | +0.27 dB [+0.14, +0.40] | +0.21 dB [+0.07, +0.26] |

## What drives Delta: placement and helmet fit (dense OPM vs Neuromag combined, intrinsic + brain)

Each child placement is compared with the adult at the same rule (the adult's counterfactual helmet has factor 1).

| child anatomy | top (primary) | centred | x-centred | top-18mm | back | counterfactual | counterfactual, x-centred |
|---|---|---|---|---|---|---|---|
| school-age size (scaled adult) | +0.44 dB [+0.34, +0.56] | +1.14 dB [+1.00, +1.25] | +0.44 dB [+0.33, +0.55] | +0.44 dB [+0.33, +0.56] | +0.85 dB [+0.65, +1.05] | -0.22 dB [-0.25, -0.18] | -0.21 dB [-0.25, -0.17] |
| 2-year size (scaled adult) | +1.07 dB [+0.91, +1.18] | +1.68 dB [+1.50, +1.82] | +1.06 dB [+0.90, +1.19] | +1.09 dB [+0.96, +1.22] | +1.21 dB [+0.91, +1.47] | -0.37 dB [-0.43, -0.33] | -0.38 dB [-0.43, -0.32] |
| 2-year template | +0.73 dB [+0.46, +0.98] | +1.83 dB [+1.46, +2.32] | +0.71 dB [+0.52, +0.86] | +0.77 dB [+0.45, +0.90] | +1.39 dB [+1.05, +1.67] | +0.24 dB [+0.13, +0.56] | -0.30 dB [-0.42, +0.02] |
| 18-month template | +0.88 dB [+0.54, +1.05] | +2.08 dB [+1.71, +2.22] | +0.75 dB [+0.48, +1.08] | +0.85 dB [+0.57, +1.24] | +1.49 dB [+1.17, +1.76] | +0.18 dB [+0.11, +0.47] | -0.16 dB [-0.39, +0.03] |
| 12-month template | +0.96 dB [+0.61, +1.51] | +2.50 dB [+2.20, +2.97] | +0.91 dB [+0.60, +1.45] | +0.85 dB [+0.54, +1.46] | +1.68 dB [+1.03, +2.26] | -0.02 dB [-0.16, +0.11] | -0.18 dB [-0.37, +0.06] |
| child A (7.8 y) | +0.30 dB [+0.11, +0.38] | +0.76 dB [+0.60, +1.18] | +0.29 dB [+0.14, +0.39] | +0.31 dB [+0.15, +0.39] | +0.44 dB [+0.28, +0.79] | -0.74 dB [-0.86, -0.60] | -0.76 dB [-0.84, -0.68] |
| child B (8.3 y) | +0.41 dB [+0.28, +0.54] | +1.57 dB [+1.16, +1.72] | +0.36 dB [+0.19, +0.45] | +0.40 dB [+0.30, +0.54] | +0.94 dB [+0.66, +1.32] | -0.23 dB [-0.49, -0.05] | -0.59 dB [-0.74, -0.42] |
| child C (8.7 y) | +0.17 dB [-0.07, +0.34] | +1.00 dB [+0.74, +1.12] | +0.09 dB [-0.03, +0.37] | +0.17 dB [-0.00, +0.36] | +0.87 dB [+0.70, +1.13] | +0.56 dB [+0.37, +0.78] | +0.04 dB [-0.20, +0.20] |

Counterfactual Delta by comparator (helmet scaled with the head; the dependence on the comparator points to the SQUID side of the change):

| child anatomy | comparator | counterfactual | counterfactual, x-centred |
|---|---|---|---|
| school-age size (scaled adult) | Neuromag combined | -0.22 dB [-0.25, -0.18] | -0.21 dB [-0.25, -0.17] |
| school-age size (scaled adult) | Neuromag grad | -0.37 dB [-0.42, -0.32] | -0.37 dB [-0.43, -0.31] |
| school-age size (scaled adult) | Neuromag mag | -0.16 dB [-0.20, -0.13] | -0.16 dB [-0.19, -0.12] |
| 2-year size (scaled adult) | Neuromag combined | -0.37 dB [-0.43, -0.33] | -0.38 dB [-0.43, -0.32] |
| 2-year size (scaled adult) | Neuromag grad | -0.60 dB [-0.68, -0.54] | -0.60 dB [-0.68, -0.53] |
| 2-year size (scaled adult) | Neuromag mag | -0.28 dB [-0.32, -0.24] | -0.28 dB [-0.33, -0.24] |
| 2-year template | Neuromag combined | +0.24 dB [+0.13, +0.56] | -0.30 dB [-0.42, +0.02] |
| 2-year template | Neuromag grad | +0.32 dB [+0.01, +0.62] | -0.42 dB [-0.68, -0.28] |
| 2-year template | Neuromag mag | +0.38 dB [+0.27, +0.46] | -0.14 dB [-0.30, +0.13] |
| 18-month template | Neuromag combined | +0.18 dB [+0.11, +0.47] | -0.16 dB [-0.39, +0.03] |
| 18-month template | Neuromag grad | +0.30 dB [+0.03, +0.55] | -0.29 dB [-0.54, -0.15] |
| 18-month template | Neuromag mag | +0.32 dB [+0.17, +0.48] | -0.07 dB [-0.28, +0.15] |
| 12-month template | Neuromag combined | -0.02 dB [-0.16, +0.11] | -0.18 dB [-0.37, +0.06] |
| 12-month template | Neuromag grad | -0.26 dB [-0.45, +0.17] | -0.36 dB [-0.67, -0.00] |
| 12-month template | Neuromag mag | +0.03 dB [-0.05, +0.18] | -0.12 dB [-0.20, +0.21] |
| child A (7.8 y) | Neuromag combined | -0.74 dB [-0.86, -0.60] | -0.76 dB [-0.84, -0.68] |
| child A (7.8 y) | Neuromag grad | -1.16 dB [-1.21, -0.96] | -1.14 dB [-1.25, -1.01] |
| child A (7.8 y) | Neuromag mag | -0.57 dB [-0.61, -0.36] | -0.58 dB [-0.64, -0.37] |
| child B (8.3 y) | Neuromag combined | -0.23 dB [-0.49, -0.05] | -0.59 dB [-0.74, -0.42] |
| child B (8.3 y) | Neuromag grad | -0.47 dB [-0.65, -0.27] | -0.94 dB [-1.04, -0.69] |
| child B (8.3 y) | Neuromag mag | -0.09 dB [-0.17, +0.10] | -0.32 dB [-0.49, -0.21] |
| child C (8.7 y) | Neuromag combined | +0.56 dB [+0.37, +0.78] | +0.04 dB [-0.20, +0.20] |
| child C (8.7 y) | Neuromag grad | +0.73 dB [+0.48, +0.87] | -0.16 dB [-0.38, +0.20] |
| child C (8.7 y) | Neuromag mag | +0.61 dB [+0.24, +0.78] | +0.07 dB [-0.11, +0.24] |

## Absolute detectability (median 20 log10 d of a 10-nAm dipole, intrinsic + brain, primary placement)

| anatomy | OPM dense | OPM matched | Neuromag combined | Neuromag grad | Neuromag mag |
|---|---|---|---|---|---|
| adult | -1.60 | -2.71 | -2.67 | -3.62 | -2.99 |
| school-age size (scaled adult) | -0.61 | -1.75 | -2.15 | -3.35 | -2.44 |
| 2-year size (scaled adult) | -0.20 | -1.21 | -2.27 | -3.65 | -2.52 |
| 2-year template | +0.92 | -0.41 | -1.19 | -2.70 | -1.52 |
| 18-month template | +1.07 | -0.30 | -1.03 | -2.66 | -1.34 |
| 12-month template | +1.86 | +0.44 | -0.62 | -2.23 | -0.97 |
| child A (7.8 y) | -1.55 | -2.86 | -2.91 | -3.92 | -3.25 |
| child B (8.3 y) | -1.41 | -2.87 | -2.70 | -3.82 | -3.04 |
| child C (8.7 y) | -2.11 | -3.25 | -3.26 | -4.36 | -3.61 |

Vertex-wise change from the adult (scaled controls; same vertex): both systems gain, the OPM more.

| child anatomy | OPM dense | Neuromag combined | Neuromag grad | Neuromag mag |
|---|---|---|---|---|
| school-age size (scaled adult) | +1.10 | +0.65 | +0.45 | +0.69 |
| 2-year size (scaled adult) | +1.67 | +0.54 | +0.10 | +0.64 |

## Extended sources (geodesic patches, primary placement, intrinsic + brain)

Median D (dense OPM vs Neuromag combined) over the patch centres (cortical, area-weighted) and the focal D at the same centres; the 20-mm patches use every third of the 300 centres. D does not depend on the moment convention.

| anatomy | radius [mm] | centres | median area [cm2] | D patch | D focal, same centres | Delta patch |
|---|---|---|---|---|---|---|
| adult | 5 | 300 | 0.69 | +1.03 | +1.03 |  |
| adult | 10 | 300 | 2.82 | +1.05 | +1.03 |  |
| adult | 20 | 100 | 11.46 | +0.83 | +1.05 |  |
| school-age size (scaled adult) | 5 | 300 | 0.70 | +1.62 | +1.53 | +0.58 |
| school-age size (scaled adult) | 10 | 300 | 2.80 | +1.59 | +1.53 | +0.54 |
| school-age size (scaled adult) | 20 | 100 | 11.39 | +1.27 | +1.66 | +0.44 |
| 2-year size (scaled adult) | 5 | 300 | 0.69 | +2.33 | +2.39 | +1.30 |
| 2-year size (scaled adult) | 10 | 300 | 2.80 | +2.26 | +2.39 | +1.22 |
| 2-year size (scaled adult) | 20 | 100 | 11.31 | +1.82 | +2.52 | +0.99 |
| 2-year template | 5 | 300 | 0.69 | +1.59 | +1.65 | +0.56 |
| 2-year template | 10 | 300 | 2.81 | +1.63 | +1.65 | +0.58 |
| 2-year template | 20 | 100 | 11.59 | +1.55 | +1.65 | +0.72 |
| 18-month template | 5 | 300 | 0.69 | +2.09 | +2.00 | +1.05 |
| 18-month template | 10 | 300 | 2.81 | +1.85 | +2.00 | +0.81 |
| 18-month template | 20 | 100 | 11.66 | +1.72 | +2.25 | +0.89 |
| 12-month template | 5 | 300 | 0.69 | +2.08 | +2.07 | +1.05 |
| 12-month template | 10 | 300 | 2.82 | +2.07 | +2.07 | +1.02 |
| 12-month template | 20 | 100 | 11.47 | +1.74 | +2.38 | +0.91 |
| child A (7.8 y) | 5 | 300 | 0.69 | +1.24 | +1.25 | +0.20 |
| child A (7.8 y) | 10 | 300 | 2.81 | +1.35 | +1.25 | +0.31 |
| child A (7.8 y) | 20 | 100 | 11.42 | +1.32 | +1.65 | +0.50 |
| child B (8.3 y) | 5 | 300 | 0.69 | +1.53 | +1.47 | +0.49 |
| child B (8.3 y) | 10 | 300 | 2.81 | +1.62 | +1.47 | +0.57 |
| child B (8.3 y) | 20 | 100 | 11.46 | +1.17 | +1.25 | +0.35 |
| child C (8.7 y) | 5 | 300 | 0.69 | +1.09 | +1.16 | +0.06 |
| child C (8.7 y) | 10 | 300 | 2.87 | +1.11 | +1.16 | +0.06 |
| child C (8.7 y) | 20 | 100 | 11.94 | +1.18 | +0.96 | +0.35 |

Absolute detectability at a fixed current density of 0.5 nAm/mm^2 (moment = density x patch area; human neocortex 0.16-0.77 nAm/mm^2, Murakami & Okada 2015): median detectability and the share of patch centres at or above the usefulness threshold (5):

| anatomy | radius [mm] | Neuromag combined | OPM dense | usable, Neuromag / OPM |
|---|---|---|---|---|
| adult | 5 | 1.7 | 2.0 | 0% / 8% |
| adult | 10 | 4.5 | 5.3 | 43% / 54% |
| adult | 20 | 9.8 | 11.1 | 91% / 92% |
| school-age size (scaled adult) | 5 | 1.7 | 2.1 | 0% / 8% |
| school-age size (scaled adult) | 10 | 4.4 | 5.4 | 41% / 53% |
| school-age size (scaled adult) | 20 | 9.2 | 11.2 | 90% / 95% |
| 2-year size (scaled adult) | 5 | 1.6 | 2.2 | 0% / 7% |
| 2-year size (scaled adult) | 10 | 4.1 | 5.4 | 31% / 54% |
| 2-year size (scaled adult) | 20 | 8.3 | 10.8 | 86% / 94% |
| 2-year template | 5 | 2.4 | 2.8 | 8% / 27% |
| 2-year template | 10 | 7.2 | 9.1 | 66% / 76% |
| 2-year template | 20 | 15.1 | 19.2 | 99% / 100% |
| 18-month template | 5 | 2.5 | 3.2 | 7% / 27% |
| 18-month template | 10 | 7.2 | 9.5 | 71% / 80% |
| 18-month template | 20 | 13.8 | 17.4 | 97% / 100% |
| 12-month template | 5 | 2.7 | 3.6 | 6% / 32% |
| 12-month template | 10 | 8.0 | 10.5 | 77% / 85% |
| 12-month template | 20 | 16.4 | 20.2 | 100% / 100% |
| child A (7.8 y) | 5 | 1.5 | 1.8 | 0% / 6% |
| child A (7.8 y) | 10 | 4.1 | 5.0 | 40% / 51% |
| child A (7.8 y) | 20 | 10.5 | 12.2 | 85% / 93% |
| child B (8.3 y) | 5 | 1.7 | 2.0 | 2% / 9% |
| child B (8.3 y) | 10 | 4.5 | 5.4 | 43% / 53% |
| child B (8.3 y) | 20 | 9.0 | 10.5 | 86% / 91% |
| child C (8.7 y) | 5 | 1.8 | 2.1 | 1% / 10% |
| child C (8.7 y) | 10 | 4.2 | 5.0 | 36% / 50% |
| child C (8.7 y) | 20 | 8.6 | 9.4 | 88% / 92% |

## Channel count: the adult's dense array subsampled to each child's site count

| child anatomy | sites | D_child | D_adult, subsampled | Delta at equal channel count |
|---|---|---|---|---|
| school-age size (scaled adult) | 174 | +1.50 dB [+1.32, +1.75] | +0.73 dB [+0.58, +0.88] | +0.63 dB [+0.50, +0.77] |
| 2-year size (scaled adult) | 155 | +2.13 dB [+1.81, +2.48] | +0.59 dB [+0.43, +0.72] | +1.41 dB [+1.23, +1.56] |
| 2-year template | 151 | +1.84 dB [+1.56, +2.19] | +0.54 dB [+0.41, +0.69] | +1.23 dB [+1.01, +1.43] |
| 18-month template | 157 | +1.89 dB [+1.65, +2.20] | +0.61 dB [+0.48, +0.74] | +1.37 dB [+0.90, +1.45] |
| 12-month template | 144 | +2.04 dB [+1.67, +2.42] | +0.46 dB [+0.34, +0.60] | +1.43 dB [+0.99, +2.17] |
| child A (7.8 y) | 155 | +1.34 dB [+1.08, +1.61] | +0.59 dB [+0.45, +0.71] | +0.66 dB [+0.43, +0.94] |
| child B (8.3 y) | 153 | +1.51 dB [+1.25, +1.78] | +0.56 dB [+0.43, +0.72] | +0.78 dB [+0.66, +0.98] |
| child C (8.7 y) | 167 | +1.19 dB [+1.00, +1.40] | +0.68 dB [+0.53, +0.83] | +0.43 dB [+0.26, +0.75] |

## Delta by depth stratum (dense OPM vs Neuromag combined, intrinsic + brain)

| child anatomy | depth [mm] | n child / adult | D_child | D_adult | Delta [95 % CI] |
|---|---|---|---|---|---|
| school-age size (scaled adult) | 0-10 | 0 / 0 | sparse | | |
| school-age size (scaled adult) | 10-15 | 450 / 222 | +3.98 | +3.41 | +0.57 [+0.11, +0.90] |
| school-age size (scaled adult) | 15-20 | 1604 / 1221 | +2.86 | +2.27 | +0.60 [+0.40, +0.86] |
| school-age size (scaled adult) | 20-25 | 1604 / 1572 | +1.77 | +1.48 | +0.29 [+0.13, +0.50] |
| school-age size (scaled adult) | 25-30 | 1144 / 1289 | +1.03 | +0.87 | +0.16 [+0.02, +0.33] |
| school-age size (scaled adult) | 30-40 | 1399 / 1539 | +0.59 | +0.45 | +0.15 [+0.03, +0.26] |
| school-age size (scaled adult) | 40-50 | 633 / 913 | +0.41 | +0.30 | +0.11 [-0.04, +0.23] |
| school-age size (scaled adult) | 50-60 | 115 / 321 | +0.74 | +0.34 | +0.40 [-0.50, +0.54] |
| school-age size (scaled adult) | 60-90 | 0 / 26 | sparse | | |
| 2-year size (scaled adult) | 0-10 | 2 / 0 | sparse | | |
| 2-year size (scaled adult) | 10-15 | 610 / 222 | +4.92 | +3.41 | +1.52 [+1.09, +1.84] |
| 2-year size (scaled adult) | 15-20 | 1768 / 1221 | +3.47 | +2.27 | +1.20 [+1.00, +1.46] |
| 2-year size (scaled adult) | 20-25 | 1610 / 1572 | +2.17 | +1.48 | +0.68 [+0.52, +0.89] |
| 2-year size (scaled adult) | 25-30 | 1055 / 1289 | +1.27 | +0.87 | +0.40 [+0.21, +0.58] |
| 2-year size (scaled adult) | 30-40 | 1301 / 1539 | +0.77 | +0.45 | +0.32 [+0.20, +0.49] |
| 2-year size (scaled adult) | 40-50 | 463 / 913 | +0.63 | +0.30 | +0.33 [+0.15, +0.51] |
| 2-year size (scaled adult) | 50-60 | 46 / 321 | +1.16 | +0.34 | +0.81 [-0.45, +0.98] |
| 2-year size (scaled adult) | 60-90 | 0 / 26 | sparse | | |
| 2-year template | 0-10 | 26 / 0 | sparse | | |
| 2-year template | 10-15 | 1488 / 222 | +4.08 | +3.41 | +0.67 [+0.18, +1.14] |
| 2-year template | 15-20 | 1666 / 1221 | +2.63 | +2.27 | +0.37 [+0.03, +0.75] |
| 2-year template | 20-25 | 1410 / 1572 | +1.65 | +1.48 | +0.17 [-0.05, +0.41] |
| 2-year template | 25-30 | 980 / 1289 | +1.04 | +0.87 | +0.17 [-0.01, +0.35] |
| 2-year template | 30-40 | 1175 / 1539 | +0.73 | +0.45 | +0.29 [+0.12, +0.55] |
| 2-year template | 40-50 | 502 / 913 | +0.54 | +0.30 | +0.24 [+0.10, +0.60] |
| 2-year template | 50-60 | 229 / 321 | +1.19 | +0.34 | +0.85 [+0.08, +1.17] |
| 2-year template | 60-90 | 83 / 26 | +2.19 | +0.58 | +1.62 [+1.44, +1.92] |
| 18-month template | 0-10 | 32 / 0 | sparse | | |
| 18-month template | 10-15 | 744 / 222 | +4.21 | +3.41 | +0.80 [+0.30, +1.38] |
| 18-month template | 15-20 | 1788 / 1221 | +2.94 | +2.27 | +0.68 [+0.30, +1.25] |
| 18-month template | 20-25 | 1539 / 1572 | +1.92 | +1.48 | +0.44 [+0.12, +0.87] |
| 18-month template | 25-30 | 1161 / 1289 | +1.31 | +0.87 | +0.44 [+0.21, +0.67] |
| 18-month template | 30-40 | 1290 / 1539 | +0.85 | +0.45 | +0.41 [+0.24, +0.63] |
| 18-month template | 40-50 | 523 / 913 | +0.54 | +0.30 | +0.24 [+0.11, +0.56] |
| 18-month template | 50-60 | 221 / 321 | +1.24 | +0.34 | +0.90 [+0.14, +1.17] |
| 18-month template | 60-90 | 65 / 26 | +2.10 | +0.58 | +1.52 [+1.42, +1.83] |
| 12-month template | 0-10 | 78 / 0 | sparse | | |
| 12-month template | 10-15 | 1744 / 222 | +4.55 | +3.41 | +1.15 [+0.45, +1.77] |
| 12-month template | 15-20 | 1768 / 1221 | +2.90 | +2.27 | +0.64 [+0.11, +1.03] |
| 12-month template | 20-25 | 1452 / 1572 | +1.78 | +1.48 | +0.30 [-0.03, +0.62] |
| 12-month template | 25-30 | 928 / 1289 | +1.09 | +0.87 | +0.22 [+0.04, +0.42] |
| 12-month template | 30-40 | 1029 / 1539 | +0.71 | +0.45 | +0.26 [+0.11, +0.45] |
| 12-month template | 40-50 | 390 / 913 | +0.58 | +0.30 | +0.28 [+0.11, +0.49] |
| 12-month template | 50-60 | 184 / 321 | +1.12 | +0.34 | +0.78 [+0.24, +0.95] |
| 12-month template | 60-90 | 19 / 26 | +1.42 | +0.58 | +0.84 [+0.77, +1.15] |
| child A (7.8 y) | 0-10 | 43 / 0 | sparse | | |
| child A (7.8 y) | 10-15 | 1188 / 222 | +3.18 | +3.41 | -0.23 [-0.61, +0.14] |
| child A (7.8 y) | 15-20 | 1699 / 1221 | +2.04 | +2.27 | -0.23 [-0.55, +0.11] |
| child A (7.8 y) | 20-25 | 1409 / 1572 | +1.20 | +1.48 | -0.28 [-0.49, -0.09] |
| child A (7.8 y) | 25-30 | 987 / 1289 | +0.74 | +0.87 | -0.13 [-0.25, -0.01] |
| child A (7.8 y) | 30-40 | 1089 / 1539 | +0.52 | +0.45 | +0.08 [-0.02, +0.25] |
| child A (7.8 y) | 40-50 | 461 / 913 | +0.45 | +0.30 | +0.15 [-0.01, +0.85] |
| child A (7.8 y) | 50-60 | 247 / 321 | +1.00 | +0.34 | +0.66 [+0.10, +1.25] |
| child A (7.8 y) | 60-90 | 15 / 26 | +1.22 | +0.58 | +0.65 [+0.58, +1.19] |
| child B (8.3 y) | 0-10 | 247 / 0 | sparse | | |
| child B (8.3 y) | 10-15 | 1550 / 222 | +3.18 | +3.41 | -0.23 [-0.58, +0.19] |
| child B (8.3 y) | 15-20 | 1538 / 1221 | +1.92 | +2.27 | -0.34 [-0.60, -0.07] |
| child B (8.3 y) | 20-25 | 1131 / 1572 | +1.12 | +1.48 | -0.36 [-0.55, -0.16] |
| child B (8.3 y) | 25-30 | 780 / 1289 | +0.73 | +0.87 | -0.14 [-0.31, +0.08] |
| child B (8.3 y) | 30-40 | 1042 / 1539 | +0.60 | +0.45 | +0.15 [+0.01, +0.40] |
| child B (8.3 y) | 40-50 | 429 / 913 | +0.78 | +0.30 | +0.48 [+0.19, +1.23] |
| child B (8.3 y) | 50-60 | 146 / 321 | +2.11 | +0.34 | +1.77 [+1.13, +2.07] |
| child B (8.3 y) | 60-90 | 3 / 26 | sparse | | |
| child C (8.7 y) | 0-10 | 11 / 0 | sparse | | |
| child C (8.7 y) | 10-15 | 842 / 222 | +3.49 | +3.41 | +0.08 [-0.32, +0.53] |
| child C (8.7 y) | 15-20 | 1602 / 1221 | +2.19 | +2.27 | -0.07 [-0.35, +0.18] |
| child C (8.7 y) | 20-25 | 1515 / 1572 | +1.28 | +1.48 | -0.20 [-0.45, -0.02] |
| child C (8.7 y) | 25-30 | 1155 / 1289 | +0.75 | +0.87 | -0.12 [-0.28, +0.05] |
| child C (8.7 y) | 30-40 | 1333 / 1539 | +0.46 | +0.45 | +0.01 [-0.11, +0.16] |
| child C (8.7 y) | 40-50 | 639 / 913 | +0.37 | +0.30 | +0.07 [-0.09, +0.49] |
| child C (8.7 y) | 50-60 | 280 / 321 | +0.70 | +0.34 | +0.36 [-0.03, +1.04] |
| child C (8.7 y) | 60-90 | 15 / 26 | +0.98 | +0.58 | +0.41 [+0.21, +1.20] |

Scaled controls, vertex-wise (homologous) Delta by the adult's depth:

| child anatomy | adult depth [mm] | n | Delta [95 % CI] |
|---|---|---|---|
| school-age size (scaled adult) | 0-10 | 0 | sparse |
| school-age size (scaled adult) | 10-15 | 175 | +1.05 [+0.57, +1.28] |
| school-age size (scaled adult) | 15-20 | 1134 | +0.89 [+0.74, +1.16] |
| school-age size (scaled adult) | 20-25 | 1564 | +0.76 [+0.62, +0.91] |
| school-age size (scaled adult) | 25-30 | 1284 | +0.52 [+0.43, +0.65] |
| school-age size (scaled adult) | 30-40 | 1533 | +0.26 [+0.20, +0.34] |
| school-age size (scaled adult) | 40-50 | 912 | +0.13 [+0.07, +0.19] |
| school-age size (scaled adult) | 50-60 | 321 | +0.14 [+0.05, +0.23] |
| school-age size (scaled adult) | 60-90 | 26 | +0.33 [-0.56, +0.40] |
| 2-year size (scaled adult) | 0-10 | 0 | sparse |
| 2-year size (scaled adult) | 10-15 | 152 | +2.40 [+1.78, +2.54] |
| 2-year size (scaled adult) | 15-20 | 1082 | +1.94 [+1.80, +2.11] |
| 2-year size (scaled adult) | 20-25 | 1551 | +1.57 [+1.44, +1.66] |
| 2-year size (scaled adult) | 25-30 | 1283 | +1.13 [+1.03, +1.22] |
| 2-year size (scaled adult) | 30-40 | 1532 | +0.64 [+0.58, +0.70] |
| 2-year size (scaled adult) | 40-50 | 909 | +0.35 [+0.29, +0.43] |
| 2-year size (scaled adult) | 50-60 | 320 | +0.33 [+0.21, +0.46] |
| 2-year size (scaled adult) | 60-90 | 26 | +0.71 [-0.30, +0.72] |

Templates and school-aged children: the pooled difference of the medians (D_child - D_adult over all targets) and the same with the child's targets reweighted to the adult's area share per depth stratum; then radial (0-30 deg) and tangential (60-90 deg) sources at matched depth (difference of the medians; '-': fewer than 10 targets):

| anatomy | pooled | depth-reweighted | area at 10-20 mm (adult) | median depth [mm] (adult) | radial 0-15 / 15-25 / 25-40 / 40-90 mm | tangential 0-15 / 15-25 / 25-40 / 40-90 mm |
|---|---|---|---|---|---|---|
| 2-year template | +0.85 | +0.32 | 42% (21%) | 21.8 (26.2) | +0.76 / +0.65 / +0.63 / +0.08 | +0.22 / +0.30 / +0.11 / +0.61 |
| 18-month template | +0.90 | +0.53 | 36% (21%) | 23.2 (26.2) | +1.19 / +1.32 / +0.88 / +0.11 | +0.38 / +0.49 / +0.20 / +0.52 |
| 12-month template | +1.05 | +0.29 | 46% (21%) | 20.6 (26.2) | +1.49 / +0.70 / +0.64 / +0.08 | +0.52 / +0.41 / +0.10 / +0.57 |
| child A (7.8 y) | +0.34 | -0.04 | 42% (21%) | 21.6 (26.2) | +0.03 / -0.57 / -0.04 / -0.14 | -0.41 / -0.16 / -0.05 / +0.39 |
| child B (8.3 y) | +0.52 | +0.05 | 47% (21%) | 20.0 (26.2) | -0.68 / -0.29 / +0.20 / +0.17 | -0.33 / -0.25 / -0.05 / +0.70 |
| child C (8.7 y) | +0.19 | -0.05 | 35% (21%) | 23.4 (26.2) | -0.03 / -0.12 / +0.01 / -0.05 | -0.13 / -0.04 / -0.06 / +0.17 |

## Delta by orientation stratum (0 deg = radial to the inner skull; dense OPM vs Neuromag combined)

| child anatomy | orientation [deg] | n child / adult | D_child | D_adult | Delta [95 % CI] |
|---|---|---|---|---|---|
| school-age size (scaled adult) | 0-30 | 676 / 707 | +1.25 | +0.84 | +0.41 [+0.20, +0.63] |
| school-age size (scaled adult) | 30-60 | 2297 / 2346 | +1.27 | +0.83 | +0.44 [+0.15, +0.73] |
| school-age size (scaled adult) | 60-90.1 | 3976 / 4050 | +1.70 | +1.13 | +0.58 [+0.22, +0.91] |
| 2-year size (scaled adult) | 0-30 | 665 / 707 | +1.77 | +0.84 | +0.93 [+0.70, +1.16] |
| 2-year size (scaled adult) | 30-60 | 2269 / 2346 | +1.84 | +0.83 | +1.02 [+0.66, +1.38] |
| 2-year size (scaled adult) | 60-90.1 | 3921 / 4050 | +2.42 | +1.13 | +1.29 [+0.86, +1.64] |
| 2-year template | 0-30 | 1219 / 707 | +2.78 | +0.84 | +1.94 [+1.33, +2.50] |
| 2-year template | 30-60 | 2293 / 2346 | +1.76 | +0.83 | +0.93 [+0.57, +1.29] |
| 2-year template | 60-90.1 | 4047 / 4050 | +1.72 | +1.13 | +0.59 [+0.25, +0.93] |
| 18-month template | 0-30 | 1026 / 707 | +2.57 | +0.84 | +1.73 [+1.17, +2.19] |
| 18-month template | 30-60 | 2398 / 2346 | +1.91 | +0.83 | +1.08 [+0.69, +1.46] |
| 18-month template | 60-90.1 | 3939 / 4050 | +1.77 | +1.13 | +0.65 [+0.34, +1.01] |
| 12-month template | 0-30 | 1294 / 707 | +3.16 | +0.84 | +2.32 [+1.46, +3.37] |
| 12-month template | 30-60 | 2461 / 2346 | +1.90 | +0.83 | +1.07 [+0.59, +1.53] |
| 12-month template | 60-90.1 | 3837 / 4050 | +1.87 | +1.13 | +0.75 [+0.35, +1.16] |
| child A (7.8 y) | 0-30 | 571 / 707 | +0.88 | +0.84 | +0.04 [-0.27, +0.47] |
| child A (7.8 y) | 30-60 | 2001 / 2346 | +1.18 | +0.83 | +0.35 [+0.03, +0.68] |
| child A (7.8 y) | 60-90.1 | 4566 / 4050 | +1.47 | +1.13 | +0.34 [-0.01, +0.65] |
| child B (8.3 y) | 0-30 | 647 / 707 | +1.30 | +0.84 | +0.46 [-0.03, +1.00] |
| child B (8.3 y) | 30-60 | 1989 / 2346 | +1.37 | +0.83 | +0.54 [+0.21, +0.90] |
| child B (8.3 y) | 60-90.1 | 4230 / 4050 | +1.62 | +1.13 | +0.50 [+0.15, +0.81] |
| child C (8.7 y) | 0-30 | 830 / 707 | +1.01 | +0.84 | +0.17 [-0.22, +0.52] |
| child C (8.7 y) | 30-60 | 2301 / 2346 | +1.00 | +0.83 | +0.18 [-0.11, +0.41] |
| child C (8.7 y) | 60-90.1 | 4261 / 4050 | +1.33 | +1.13 | +0.21 [-0.08, +0.51] |

## Placement, counterfactual helmet and sensitivity (median D, dense OPM vs Neuromag combined, intrinsic + brain)

| anatomy | centred | top | back | x+5mm | x-5mm | y+5mm | y-5mm | pitch+10deg | pitch-10deg | roll+5deg | roll-5deg | yaw+10deg | yaw-10deg | x-centred | top-18mm | counterfactual | counterfactual_x-centred |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| adult | +1.27 | +1.00 | +1.13 | +1.16 | +1.27 | +1.28 | +0.98 | +1.02 | +1.43 | +1.26 | +1.19 | +1.03 | +1.22 | +1.01 | +0.87 | +1.27 | +1.28 |
| school-age size (scaled adult) | +2.52 | +1.50 | +2.02 | +1.49 | +1.52 | +1.58 | +1.47 | +1.50 | +2.56 | +1.46 | +1.71 | +1.54 | +1.50 | +1.51 | +1.38 | +0.97 | +0.97 |
| 2-year size (scaled adult) | +3.05 | +2.13 | +2.34 | +2.00 | +2.15 | +2.23 | +2.05 | +1.69 | +3.04 | +1.93 | +2.50 | +2.26 | +2.20 | +2.13 | +2.03 | +0.73 | +0.73 |
| 2-year template | +3.08 | +1.84 | +2.51 | +2.07 | +1.75 | +1.84 | +1.85 | +1.78 | +2.00 | +1.79 | +2.04 | +1.90 | +1.78 | +1.76 | +1.70 | +1.55 | +1.00 |
| 18-month template | +3.13 | +1.89 | +2.48 | +2.05 | +1.80 | +1.86 | +1.82 | +1.80 | +2.03 | +1.80 | +2.04 | +1.90 | +1.81 | +1.81 | +1.75 | +1.47 | +1.02 |
| 12-month template | +3.65 | +2.04 | +2.79 | +2.13 | +2.03 | +2.07 | +2.06 | +1.98 | +2.22 | +2.01 | +2.11 | +2.07 | +1.99 | +1.98 | +1.85 | +1.16 | +0.99 |
| child A (7.8 y) | +2.07 | +1.34 | +1.40 | +1.42 | +1.40 | +1.39 | +1.33 | +1.25 | +1.50 | +1.45 | +1.42 | +1.40 | +1.31 | +1.34 | +1.24 | +0.62 | +0.62 |
| child B (8.3 y) | +2.63 | +1.51 | +1.76 | +1.77 | +1.48 | +1.53 | +1.45 | +1.40 | +1.96 | +1.45 | +1.70 | +1.55 | +1.45 | +1.47 | +1.40 | +1.06 | +0.78 |
| child C (8.7 y) | +2.08 | +1.19 | +1.83 | +1.24 | +2.00 (infeasible) | +1.18 | +1.19 | +1.05 | +1.87 | +1.91 | +1.59 | +1.33 | +1.32 | +1.25 | +1.08 | +1.73 | +1.33 |

Infeasible: a magnetometer coil centre closer to the scalp than the 18-mm Dewar spacing (not a possible position; its D is listed for completeness and left out of the placement ranges).

| anatomy | opm_asd_7fT | opm_asd_10fT | opm_asd_15fT | opm_asd_20fT | opm_asd_30fT | background_x0.5 | background_x2 | bem1 |
|---|---|---|---|---|---|---|---|---|
| adult | +2.05 | +1.58 | +1.00 | +0.56 | -0.05 | +0.97 | +1.01 | +0.97 |
| school-age size (scaled adult) | +2.52 | +2.10 | +1.50 | +1.05 | +0.37 | +1.50 | +1.49 | +1.49 |
| 2-year size (scaled adult) | +3.15 | +2.71 | +2.13 | +1.68 | +0.97 | +2.15 | +2.08 | +2.10 |
| 2-year template | +2.81 | +2.40 | +1.84 | +1.38 | +0.61 | +1.86 | +1.79 | +1.80 |
| 18-month template | +2.98 | +2.53 | +1.89 | +1.39 | +0.61 | +1.86 | +1.90 | +1.81 |
| 12-month template | +3.06 | +2.64 | +2.04 | +1.54 | +0.77 | +2.05 | +1.98 | +2.03 |
| child A (7.8 y) | +1.89 | +1.66 | +1.34 | +1.06 | +0.55 | +1.43 | +1.21 | +1.33 |
| child B (8.3 y) | +2.02 | +1.82 | +1.51 | +1.23 | +0.72 | +1.64 | +1.35 | +1.45 |
| child C (8.7 y) | +1.95 | +1.65 | +1.19 | +0.81 | +0.27 | +1.22 | +1.12 | +1.16 |

Difference of these medians, child minus adult, with the same variant applied to both (a sensitivity of the medians, not the paired Delta estimator; the background variants scale the adult too):

| child anatomy | opm_asd_7fT | opm_asd_10fT | opm_asd_15fT | opm_asd_20fT | opm_asd_30fT | background_x0.5 | background_x2 | bem1 |
|---|---|---|---|---|---|---|---|---|
| school-age size (scaled adult) | +0.48 | +0.52 | +0.51 | +0.49 | +0.42 | +0.53 | +0.48 | +0.51 |
| 2-year size (scaled adult) | +1.11 | +1.13 | +1.14 | +1.12 | +1.03 | +1.18 | +1.07 | +1.13 |
| 2-year template | +0.76 | +0.82 | +0.85 | +0.82 | +0.66 | +0.89 | +0.78 | +0.83 |
| 18-month template | +0.93 | +0.94 | +0.90 | +0.83 | +0.66 | +0.89 | +0.89 | +0.84 |
| 12-month template | +1.01 | +1.05 | +1.05 | +0.98 | +0.82 | +1.08 | +0.97 | +1.06 |
| child A (7.8 y) | -0.15 | +0.08 | +0.34 | +0.50 | +0.60 | +0.46 | +0.20 | +0.35 |
| child B (8.3 y) | -0.03 | +0.24 | +0.52 | +0.67 | +0.77 | +0.67 | +0.34 | +0.48 |
| child C (8.7 y) | -0.10 | +0.06 | +0.19 | +0.25 | +0.32 | +0.25 | +0.11 | +0.19 |

## Regions that gain or lose with the placement (median D by lobe, dense OPM vs Neuromag combined, intrinsic + brain)

| anatomy | placement | frontal | parietal | temporal | occipital | cingulate | insula |
|---|---|---|---|---|---|---|---|
| adult | centred | +1.77 | +1.14 | +1.34 | +0.93 | +0.40 | +0.73 |
| adult | top | +1.40 | +0.84 | +1.17 | +0.79 | +0.25 | +0.56 |
| adult | x-centred | +1.43 | +0.87 | +1.19 | +0.80 | +0.25 | +0.52 |
| adult | back | +1.91 | +0.94 | +1.17 | +0.47 | +0.37 | +0.65 |
| adult | counterfactual | +1.77 | +1.14 | +1.34 | +0.93 | +0.40 | +0.73 |
| school-age size (scaled adult) | centred | +3.14 | +2.57 | +2.53 | +1.95 | +0.88 | +1.31 |
| school-age size (scaled adult) | top | +1.78 | +1.37 | +1.91 | +1.48 | +0.35 | +0.83 |
| school-age size (scaled adult) | x-centred | +1.79 | +1.37 | +1.90 | +1.48 | +0.35 | +0.81 |
| school-age size (scaled adult) | back | +3.43 | +1.90 | +1.99 | +0.74 | +0.74 | +1.12 |
| school-age size (scaled adult) | counterfactual | +1.40 | +0.96 | +0.94 | +0.59 | +0.28 | +0.50 |
| 2-year size (scaled adult) | centred | +3.77 | +3.06 | +3.14 | +2.39 | +1.10 | +1.62 |
| 2-year size (scaled adult) | top | +2.58 | +2.03 | +2.50 | +1.93 | +0.59 | +1.20 |
| 2-year size (scaled adult) | x-centred | +2.58 | +2.03 | +2.48 | +1.92 | +0.59 | +1.19 |
| 2-year size (scaled adult) | back | +4.02 | +2.12 | +2.32 | +0.76 | +0.91 | +1.36 |
| 2-year size (scaled adult) | counterfactual | +1.19 | +0.67 | +0.67 | +0.39 | +0.20 | +0.40 |
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
| 12-month template | counterfactual | +2.06 | +0.85 | +1.24 | +0.62 | +0.44 | +0.71 |
| child A (7.8 y) | centred | +2.42 | +1.91 | +2.25 | +2.26 | +0.54 | +0.60 |
| child A (7.8 y) | top | +1.44 | +0.94 | +1.90 | +1.91 | +0.35 | +0.41 |
| child A (7.8 y) | x-centred | +1.44 | +0.94 | +1.90 | +1.91 | +0.35 | +0.41 |
| child A (7.8 y) | back | +2.75 | +1.11 | +1.51 | +0.51 | +0.48 | +0.49 |
| child A (7.8 y) | counterfactual | +0.83 | +0.40 | +0.77 | +0.83 | +0.22 | +0.17 |
| child B (8.3 y) | centred | +2.69 | +2.84 | +2.92 | +2.79 | +0.87 | +1.07 |
| child B (8.3 y) | top | +1.34 | +1.18 | +2.33 | +2.08 | +0.46 | +0.72 |
| child B (8.3 y) | x-centred | +1.30 | +1.12 | +2.29 | +2.03 | +0.45 | +0.69 |
| child B (8.3 y) | back | +2.97 | +2.13 | +1.87 | +0.84 | +0.73 | +0.89 |
| child B (8.3 y) | counterfactual | +1.10 | +0.98 | +1.30 | +1.26 | +0.41 | +0.37 |
| child C (8.7 y) | centred | +2.50 | +2.46 | +1.73 | +2.35 | +0.62 | +0.69 |
| child C (8.7 y) | top | +1.20 | +1.22 | +1.23 | +1.71 | +0.26 | +0.29 |
| child C (8.7 y) | x-centred | +1.19 | +1.26 | +1.41 | +1.82 | +0.28 | +0.35 |
| child C (8.7 y) | back | +2.68 | +2.12 | +1.41 | +1.53 | +0.60 | +0.61 |
| child C (8.7 y) | counterfactual | +2.10 | +2.07 | +1.33 | +1.96 | +0.52 | +0.53 |

## Source-to-sensor distance (median, mm) by depth below the scalp

| anatomy | sensors | 0-10 | 10-15 | 15-20 | 20-25 | 25-30 | 30-40 | 40-50 | 50-60 | 60-90 |
|---|---|---|---|---|---|---|---|---|---|---|
| adult | squid:top | - | 45 | 48 | 52 | 57 | 64 | 74 | 82 | 87 |
| adult | opm_dense | - | 22 | 27 | 31 | 36 | 43 | 53 | 63 | 69 |
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
| child A (7.8 y) | squid:top | 46 | 47 | 52 | 57 | 62 | 69 | 76 | 83 | 85 |
| child A (7.8 y) | opm_dense | 18 | 22 | 26 | 31 | 35 | 43 | 52 | 62 | 68 |
| child B (8.3 y) | squid:top | 44 | 49 | 54 | 59 | 64 | 72 | 77 | 82 | 89 |
| child B (8.3 y) | opm_dense | 18 | 21 | 26 | 31 | 35 | 42 | 52 | 61 | 68 |
| child C (8.7 y) | squid:top | 45 | 47 | 50 | 55 | 60 | 67 | 76 | 83 | 91 |
| child C (8.7 y) | opm_dense | 18 | 22 | 26 | 31 | 36 | 43 | 52 | 62 | 68 |

## Noise composition (median per-channel RMS in the band; magnetometers and OPM in fT, gradiometers in fT/cm)

| anatomy | channels | brain background | intrinsic |
|---|---|---|---|
| adult | squid:top/grad | 41.2 | 21.3 |
| adult | squid:top/mag | 202.4 | 20.7 |
| adult | squid:centred/grad | 37.1 | 21.3 |
| adult | squid:centred/mag | 192.1 | 20.7 |
| adult | opm_dense/mag | 494.4 | 88.8 |
| adult | opm_matched/mag | 505.5 | 88.8 |
| school-age size (scaled adult) | squid:top/grad | 29.8 | 21.3 |
| school-age size (scaled adult) | squid:top/mag | 149.8 | 20.7 |
| school-age size (scaled adult) | squid:centred/grad | 21.4 | 21.3 |
| school-age size (scaled adult) | squid:centred/mag | 131.0 | 20.7 |
| school-age size (scaled adult) | opm_dense/mag | 525.0 | 88.8 |
| school-age size (scaled adult) | opm_matched/mag | 530.6 | 88.8 |
| 2-year size (scaled adult) | squid:top/grad | 23.4 | 21.3 |
| 2-year size (scaled adult) | squid:top/mag | 124.8 | 20.7 |
| 2-year size (scaled adult) | squid:centred/grad | 16.6 | 21.3 |
| 2-year size (scaled adult) | squid:centred/mag | 108.0 | 20.7 |
| 2-year size (scaled adult) | opm_dense/mag | 537.1 | 88.8 |
| 2-year size (scaled adult) | opm_matched/mag | 543.8 | 88.8 |
| 2-year template | squid:top/grad | 23.4 | 21.3 |
| 2-year template | squid:top/mag | 119.7 | 20.7 |
| 2-year template | squid:centred/grad | 16.4 | 21.3 |
| 2-year template | squid:centred/mag | 100.1 | 20.7 |
| 2-year template | opm_dense/mag | 488.5 | 88.8 |
| 2-year template | opm_matched/mag | 483.9 | 88.8 |
| 18-month template | squid:top/grad | 21.1 | 21.3 |
| 18-month template | squid:top/mag | 110.1 | 20.7 |
| 18-month template | squid:centred/grad | 15.7 | 21.3 |
| 18-month template | squid:centred/mag | 94.3 | 20.7 |
| 18-month template | opm_dense/mag | 447.6 | 88.8 |
| 18-month template | opm_matched/mag | 462.7 | 88.8 |
| 12-month template | squid:top/grad | 18.5 | 21.3 |
| 12-month template | squid:top/mag | 98.2 | 20.7 |
| 12-month template | squid:centred/grad | 14.1 | 21.3 |
| 12-month template | squid:centred/mag | 88.1 | 20.7 |
| 12-month template | opm_dense/mag | 474.2 | 88.8 |
| 12-month template | opm_matched/mag | 470.1 | 88.8 |
| child A (7.8 y) | squid:top/grad | 39.6 | 21.3 |
| child A (7.8 y) | squid:top/mag | 193.9 | 20.7 |
| child A (7.8 y) | squid:centred/grad | 30.3 | 21.3 |
| child A (7.8 y) | squid:centred/mag | 171.6 | 20.7 |
| child A (7.8 y) | opm_dense/mag | 710.5 | 88.8 |
| child A (7.8 y) | opm_matched/mag | 696.3 | 88.8 |
| child B (8.3 y) | squid:top/grad | 33.6 | 21.3 |
| child B (8.3 y) | squid:top/mag | 168.2 | 20.7 |
| child B (8.3 y) | squid:centred/grad | 23.0 | 21.3 |
| child B (8.3 y) | squid:centred/mag | 136.7 | 20.7 |
| child B (8.3 y) | opm_dense/mag | 751.3 | 88.8 |
| child B (8.3 y) | opm_matched/mag | 743.4 | 88.8 |
| child C (8.7 y) | squid:top/grad | 36.4 | 21.3 |
| child C (8.7 y) | squid:top/mag | 174.1 | 20.7 |
| child C (8.7 y) | squid:centred/grad | 25.2 | 21.3 |
| child C (8.7 y) | squid:centred/mag | 149.5 | 20.7 |
| child C (8.7 y) | opm_dense/mag | 616.6 | 88.8 |
| child C (8.7 y) | opm_matched/mag | 626.6 | 88.8 |

## Usefulness (dense OPM vs Neuromag combined, intrinsic + brain): share of usable cortical area

A source counts as usable when its detectability reaches 5 at the reference moment (an operational choice, not a clinical standard).

| anatomy | moment [nAm] | both | OPM only | SQUID only | neither |
|---|---|---|---|---|---|
| adult | 20 | 0.00 | 0.04 | 0.00 | 0.96 |
| adult | 50 | 0.35 | 0.07 | 0.00 | 0.57 |
| adult | 100 | 0.66 | 0.02 | 0.00 | 0.32 |
| adult | 200 | 0.89 | 0.01 | 0.00 | 0.10 |
| school-age size (scaled adult) | 20 | 0.00 | 0.07 | 0.00 | 0.93 |
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
| child A (7.8 y) | 20 | 0.00 | 0.06 | 0.00 | 0.94 |
| child A (7.8 y) | 50 | 0.33 | 0.10 | 0.00 | 0.57 |
| child A (7.8 y) | 100 | 0.65 | 0.04 | 0.01 | 0.30 |
| child A (7.8 y) | 200 | 0.86 | 0.02 | 0.01 | 0.11 |
| child B (8.3 y) | 20 | 0.00 | 0.08 | 0.00 | 0.92 |
| child B (8.3 y) | 50 | 0.34 | 0.10 | 0.00 | 0.56 |
| child B (8.3 y) | 100 | 0.65 | 0.04 | 0.01 | 0.29 |
| child B (8.3 y) | 200 | 0.88 | 0.02 | 0.01 | 0.09 |
| child C (8.7 y) | 20 | 0.00 | 0.06 | 0.00 | 0.94 |
| child C (8.7 y) | 50 | 0.31 | 0.09 | 0.00 | 0.59 |
| child C (8.7 y) | 100 | 0.63 | 0.04 | 0.01 | 0.32 |
| child C (8.7 y) | 200 | 0.86 | 0.02 | 0.01 | 0.11 |

## Notes

- D = 20 log10(d_OPM / d_SQUID) of a 10-nAm cortical-normal dipole (known-topography detectability with the oracle noise covariance; independent of the moment). Delta = D_child - D_adult. A positive Delta is an increase in relative OPM performance under these matching assumptions; it does not by itself mean that OPM beats SQUID in the child.
- Scaled controls: the adult's vertices, so Delta is vertex-wise. The absolute 4-mm usable-source rule drops 154 and 248 superficial adult targets in the scaled copies, so D_child and D_adult are medians over slightly different target sets while Delta uses the common vertices. The templates and the school-aged children: no vertex correspondence; Delta is computed per Desikan-Killiany parcel and per declared depth/orientation stratum from area-weighted medians.
- School-aged children (OpenNeuro ds005234, typically developing, 7.8-8.7 years): their own white surfaces, aparc labels and MRI scalp; the dataset's watershed inner skull lies just below the scalp, so the skull is modelled (A-BEM-CHILD: the inner skull moved to 8 mm below the scalp where shallower, at least 2 mm from the cortex; the outer skull halfway to the scalp), and the fiducials are the adult's transferred by a cortex-to-cortex similarity fit (A-G3-FID). Individual children, not a population; no cortex maps are drawn for them (their inflated surfaces were not obtained).
- Intervals: bootstrap over parcels of one anatomy (or of each anatomy, for between-anatomy strata); they do not include between-subject variability. Three average templates of one database (2-year template, 18-month template, 12-month template) are not a population: template results are conditional simulations.
- Every child array uses the adult's conventions: background moment variance per unit cortical area, room field, intrinsic noise, sensor sizes and the 3-layer BEM conductivities; only geometry changes. Both systems' detectability rises in the smaller heads, the OPM's more (absolute detectability table), by different routes: the on-scalp OPM sees more signal from a cortex that is closer in absolute terms at about the same brain noise, while the SQUIDs' brain noise falls (the cortex is farther from the fixed helmet and, with the background fixed per unit area, smaller) more than their signal. The templates' averaged white surfaces are smoother than an individual cortex (usable area 1,062, 975, 895 cm^2 for the 2-year template, 18-month template, 12-month template vs 1,878 cm^2 for the adult), which lowers their background power and their patch cancellation further; scaling the background variance x0.5 or x2 leaves D_child almost unchanged.
- Placements are chosen from the scalp and helmet geometry only. Under the adult's measured pose a head with other fiducials need not be centred laterally; 'x-centred' shifts each head along device x to equal left/right median gaps before the top contact (shift: adult -1.5 mm, school-age size (scaled adult) -0.5 mm, 2-year size (scaled adult) -0.5 mm, 2-year template -6.5 mm, 18-month template -5.5 mm, 12-month template -2.5 mm, child A (7.8 y) +0.0 mm, child B (8.3 y) -2.5 mm, child C (8.7 y) +3.5 mm; negative = to the left), and 'counterfactual_x-centred' scales the helmet about that laterally centred head. The counterfactual helmet (scaled with the head) is a mechanistic control, not a pediatric SQUID system.
- Targets on the medial wall (FreeSurfer 'unknown': the cut through the corpus callosum and midbrain, not cortex) are left out of every summary; they would otherwise dominate the deepest strata.
