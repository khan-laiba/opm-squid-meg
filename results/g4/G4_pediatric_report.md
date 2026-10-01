# G4 pediatric: IED detection and bounded localization in the fixed adult helmet (NEW)

Same framework, configuration and seeds as the adult (configs/g4_epilepsy.toml); the child anatomy, arrays and placement come from G3B (Neuromag at top contact, refitted OPM arrays). Thresholds are calibrated on each anatomy's own null data. p-values are uncorrected; the location is the statistical unit.

## Strength for 50 % detection [nAm] (focal; practical detector at 1 false event/min)

| anatomy | detector | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|---|
| adult | squid/combined | 53 [44-66] | 89 [67-122] | 132 [113-184] | 290 [open] |
| adult | squid/grad | 56 [48-67] | 100 [72-138] | 140 [120-202] | none [open] |
| adult | squid/mag | 51 [40-65] | 89 [68-121] | 132 [112-177] | 296 [open] |
| adult | opm_matched/opm | 41 [32-59] | 88 [64-136] | 137 [114-193] | none [open] |
| adult | opm_dense/opm | 35 [28-50] | 82 [64-106] | 127 [105-172] | 271 [open] |
| school | squid/combined | 52 [40-68] | 90 [69-113] | 150 [117-232] | 275 [open] |
| school | squid/grad | 54 [43-66] | 101 [73-135] | 182 [139-269] | none [open] |
| school | squid/mag | 55 [43-67] | 90 [68-115] | 150 [114-226] | 271 [open] |
| school | opm_matched/opm | 46 [34-66] | 96 [72-114] | 155 [122-258] | 271 [open] |
| school | opm_dense/opm | 37 [29-51] | 85 [65-107] | 151 [119-219] | 245 [open] |
| infant2yr | squid/combined | 45 [35-54] | 71 [61-96] | 140 [118-185] | 243 [207-292] |
| infant2yr | squid/grad | 51 [44-57] | 80 [64-107] | 160 [129-218] | 279 [open] |
| infant2yr | squid/mag | 45 [35-53] | 78 [63-103] | 139 [116-187] | 240 [203-289] |
| infant2yr | opm_matched/opm | 35 [30-45] | 70 [58-92] | 150 [117-202] | 240 [208-295] |
| infant2yr | opm_dense/opm | 30 [27-35] | 66 [58-88] | 138 [110-182] | 229 [196-279] |

## Paired OPM dense vs Neuromag combined (practical detector, 1 false event/min; locations favouring OPM / SQUID, sign-flip p, S50 ratio SQUID/OPM [95 % CI])

| anatomy | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|
| adult | 16/0, p 3.05e-05, 1.51 [1.25-1.66] | 6/2, p 0.227, 1.08 [1.00-1.20] | 6/1, p 0.281, 1.04 [0.99-1.12] | 6/1, p 0.125, 1.07 [1.00-1.16] |
| school | 16/0, p 3.05e-05, 1.42 [1.20-1.57] | 6/3, p 0.285, 1.06 [0.96-1.21] | 4/3, p 0.766, 1.00 [0.93-1.10] | 6/1, p 0.109, 1.12 [1.04-1.22] |
| infant2yr | 16/2, p 0.000465, 1.50 [1.26-1.68] | 7/1, p 0.148, 1.07 [0.97-1.23] | 6/3, p 0.617, 1.01 [0.92-1.12] | 7/3, p 0.217, 1.06 [0.97-1.15] |

## Localization (median error [mm]; joint detection + localization within 10 mm)

| anatomy | array | condition | dSPM error (all) | ECD error (detected) | detected | joint dSPM | joint ECD |
|---|---|---|---|---|---|---|---|
| adult | squid | focal 80nAm | 53.0 | 5.1 | 0.29 | 0.00 | 0.25 |
| adult | squid | focal 320nAm | 14.0 | 4.3 | 0.92 | 0.21 | 0.71 |
| adult | squid | patch 80nAm | 75.2 | 16.6 | 0.12 | 0.00 | 0.00 |
| adult | squid | patch 320nAm | 19.7 | 7.2 | 0.83 | 0.17 | 0.58 |
| adult | opm_matched | focal 80nAm | 33.8 | 6.7 | 0.38 | 0.17 | 0.29 |
| adult | opm_matched | focal 320nAm | 12.4 | 4.8 | 0.92 | 0.38 | 0.62 |
| adult | opm_matched | patch 80nAm | 53.4 | 30.6 | 0.08 | 0.04 | 0.04 |
| adult | opm_matched | patch 320nAm | 15.1 | 8.5 | 0.83 | 0.38 | 0.58 |
| adult | opm_dense | focal 80nAm | 25.9 | 3.7 | 0.38 | 0.17 | 0.33 |
| adult | opm_dense | focal 320nAm | 11.7 | 4.2 | 0.92 | 0.33 | 0.75 |
| adult | opm_dense | patch 80nAm | 61.1 | 17.4 | 0.17 | 0.04 | 0.00 |
| adult | opm_dense | patch 320nAm | 12.4 | 6.7 | 0.83 | 0.42 | 0.58 |
| school | squid | focal 80nAm | 33.3 | 5.6 | 0.42 | 0.04 | 0.38 |
| school | squid | focal 320nAm | 15.9 | 7.8 | 0.92 | 0.25 | 0.58 |
| school | squid | patch 80nAm | 63.4 | - | 0.00 | 0.00 | 0.00 |
| school | squid | patch 320nAm | 15.7 | 8.8 | 0.75 | 0.25 | 0.50 |
| school | opm_matched | focal 80nAm | 19.6 | 4.7 | 0.38 | 0.08 | 0.38 |
| school | opm_matched | focal 320nAm | 10.5 | 5.2 | 0.92 | 0.38 | 0.67 |
| school | opm_matched | patch 80nAm | 46.6 | 11.2 | 0.04 | 0.04 | 0.00 |
| school | opm_matched | patch 320nAm | 13.2 | 9.1 | 0.83 | 0.33 | 0.42 |
| school | opm_dense | focal 80nAm | 21.9 | 4.9 | 0.46 | 0.04 | 0.38 |
| school | opm_dense | focal 320nAm | 12.0 | 6.3 | 0.92 | 0.29 | 0.62 |
| school | opm_dense | patch 80nAm | 48.2 | - | 0.00 | 0.00 | 0.00 |
| school | opm_dense | patch 320nAm | 14.6 | 9.3 | 0.83 | 0.33 | 0.54 |
| infant2yr | squid | focal 80nAm | 19.0 | 6.8 | 0.42 | 0.12 | 0.33 |
| infant2yr | squid | focal 320nAm | 14.8 | 4.5 | 0.83 | 0.29 | 0.67 |
| infant2yr | squid | patch 80nAm | 38.0 | 6.7 | 0.33 | 0.04 | 0.25 |
| infant2yr | squid | patch 320nAm | 13.6 | 6.4 | 0.88 | 0.42 | 0.75 |
| infant2yr | opm_matched | focal 80nAm | 22.3 | 6.5 | 0.50 | 0.33 | 0.42 |
| infant2yr | opm_matched | focal 320nAm | 8.9 | 4.7 | 0.88 | 0.54 | 0.75 |
| infant2yr | opm_matched | patch 80nAm | 39.9 | 4.5 | 0.29 | 0.12 | 0.25 |
| infant2yr | opm_matched | patch 320nAm | 10.0 | 5.0 | 0.96 | 0.50 | 0.67 |
| infant2yr | opm_dense | focal 80nAm | 16.6 | 5.3 | 0.50 | 0.29 | 0.38 |
| infant2yr | opm_dense | focal 320nAm | 7.5 | 5.8 | 0.92 | 0.54 | 0.54 |
| infant2yr | opm_dense | patch 80nAm | 39.0 | 6.2 | 0.33 | 0.12 | 0.25 |
| infant2yr | opm_dense | patch 320nAm | 10.1 | 7.2 | 1.00 | 0.50 | 0.71 |

Simulated IED-source recovery does not identify an epileptogenic zone or establish surgical benefit. One template is not a population; the scaled adult is a size-only control.
