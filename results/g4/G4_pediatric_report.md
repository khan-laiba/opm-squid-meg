# G4 pediatric: IED detection and bounded localization in the fixed adult helmet (NEW)

Same framework, configuration and seeds as the adult (configs/g4_epilepsy.toml); the child anatomy, arrays and placement come from G3B (Neuromag at top contact, refitted OPM arrays). The adult rows are the frozen adult study (Neuromag at its measured head position; top contact would raise it by 5.5 mm). Thresholds are calibrated on each anatomy's own null data. p-values are uncorrected (24 paired detection comparisons per OPM array and anatomy); the location is the statistical unit.

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
| size2yr | squid/combined | 59 [49-77] | 106 [78-135] | 202 [153-264] | 273 [open] |
| size2yr | squid/grad | 65 [53-89] | 110 [87-135] | 217 [172-287] | none [open] |
| size2yr | squid/mag | 59 [48-78] | 107 [80-134] | 207 [156-273] | 273 [open] |
| size2yr | opm_matched/opm | 38 [31-62] | 96 [71-116] | 187 [144-239] | 277 [open] |
| size2yr | opm_dense/opm | 36 [29-50] | 94 [68-126] | 171 [135-229] | 243 [open] |
| infant2yr | squid/combined | 45 [35-54] | 71 [61-96] | 140 [118-185] | 243 [207-292] |
| infant2yr | squid/grad | 51 [44-57] | 80 [64-107] | 160 [129-218] | 279 [open] |
| infant2yr | squid/mag | 45 [35-53] | 78 [63-103] | 139 [116-187] | 240 [203-289] |
| infant2yr | opm_matched/opm | 35 [30-45] | 70 [58-92] | 150 [117-202] | 240 [208-295] |
| infant2yr | opm_dense/opm | 30 [27-35] | 66 [58-88] | 138 [110-182] | 229 [196-279] |
| infant18mo | squid/combined | 36 [31-56] | 62 [54-76] | 124 [107-155] | none [open] |
| infant18mo | squid/grad | 47 [34-67] | 83 [61-118] | 144 [122-189] | none [open] |
| infant18mo | squid/mag | 37 [31-55] | 64 [54-95] | 125 [110-160] | none [open] |
| infant18mo | opm_matched/opm | 33 [26-48] | 60 [42-95] | 121 [106-148] | none [open] |
| infant18mo | opm_dense/opm | 27 [21-40] | 58 [43-88] | 118 [99-152] | 282 [open] |
| infant12mo | squid/combined | 40 [33-56] | 64 [52-91] | 123 [100-165] | 279 [open] |
| infant12mo | squid/grad | 49 [39-59] | 70 [58-101] | 143 [116-196] | 320 [open] |
| infant12mo | squid/mag | 39 [32-54] | 65 [53-93] | 119 [99-149] | 264 [open] |
| infant12mo | opm_matched/opm | 32 [28-46] | 62 [51-84] | 119 [95-160] | 266 [open] |
| infant12mo | opm_dense/opm | 29 [21-40] | 62 [50-87] | 115 [88-160] | 259 [open] |

## Held-out false events per minute at the 1-per-minute thresholds, and sensitivity at a matched held-out rate

Sensitivity for 40-nAm spikes at 10-30 mm with each detector's threshold set on the held-out null to 1 false event per minute (the frozen thresholds give 0.5-1.7 per minute, unequal between arrays).

| anatomy | detector | held-out rate at the frozen threshold | sensitivity at a matched 1 per minute |
|---|---|---|---|
| adult | squid/combined | 0.95 | 0.17 |
| adult | squid/grad | 0.70 | 0.15 |
| adult | squid/mag | 1.10 | 0.17 |
| adult | opm_matched/opm | 0.50 | 0.29 |
| adult | opm_dense/opm | 0.90 | 0.31 |
| school | squid/combined | 0.55 | 0.18 |
| school | squid/grad | 1.55 | 0.14 |
| school | squid/mag | 0.85 | 0.13 |
| school | opm_matched/opm | 0.80 | 0.23 |
| school | opm_dense/opm | 1.00 | 0.31 |
| size2yr | squid/combined | 0.40 | 0.13 |
| size2yr | squid/grad | 1.10 | 0.10 |
| size2yr | squid/mag | 0.55 | 0.13 |
| size2yr | opm_matched/opm | 0.60 | 0.29 |
| size2yr | opm_dense/opm | 0.90 | 0.31 |
| infant2yr | squid/combined | 0.65 | 0.23 |
| infant2yr | squid/grad | 0.60 | 0.18 |
| infant2yr | squid/mag | 1.70 | 0.20 |
| infant2yr | opm_matched/opm | 1.55 | 0.34 |
| infant2yr | opm_dense/opm | 1.25 | 0.43 |
| infant18mo | squid/combined | 0.95 | 0.38 |
| infant18mo | squid/grad | 0.75 | 0.31 |
| infant18mo | squid/mag | 1.10 | 0.35 |
| infant18mo | opm_matched/opm | 1.45 | 0.44 |
| infant18mo | opm_dense/opm | 1.05 | 0.50 |
| infant12mo | squid/combined | 1.25 | 0.32 |
| infant12mo | squid/grad | 1.15 | 0.23 |
| infant12mo | squid/mag | 1.35 | 0.34 |
| infant12mo | opm_matched/opm | 1.05 | 0.42 |
| infant12mo | opm_dense/opm | 0.70 | 0.46 |

## Paired OPM dense vs Neuromag combined (practical detector, 1 false event/min; locations favouring OPM / SQUID, sign-flip p, S50 ratio SQUID/OPM [95 % CI])

| anatomy | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|
| adult | 16/0, p 3.05e-05, 1.51 [1.25-1.66] | 6/2, p 0.227, 1.08 [1.00-1.20] | 6/1, p 0.281, 1.04 [0.99-1.12] | 6/1, p 0.125, 1.07 [1.00-1.16] |
| school | 16/0, p 3.05e-05, 1.42 [1.20-1.57] | 6/3, p 0.285, 1.06 [0.96-1.21] | 4/3, p 0.766, 1.00 [0.93-1.10] | 6/1, p 0.109, 1.12 [1.04-1.22] |
| size2yr | 16/0, p 3.05e-05, 1.64 [1.40-1.85] | 7/3, p 0.189, 1.12 [0.95-1.37] | 10/1, p 0.00879, 1.18 [1.06-1.31] | 8/1, p 0.0312, 1.12 [1.03-1.24] |
| infant2yr | 16/2, p 0.000465, 1.50 [1.26-1.68] | 7/1, p 0.148, 1.07 [0.97-1.23] | 6/3, p 0.617, 1.01 [0.92-1.12] | 7/3, p 0.217, 1.06 [0.97-1.15] |
| infant18mo | 16/0, p 3.05e-05, 1.33 [1.16-1.62] | 6/1, p 0.0938, 1.06 [0.90-1.26] | 8/3, p 0.183, 1.05 [0.98-1.13] | 7/1, p 0.0469, n/a [1.03-1.29] |
| infant12mo | 14/1, p 0.000427, 1.39 [1.19-1.68] | 5/3, p 0.562, 1.04 [0.96-1.17] | 7/4, p 0.432, 1.07 [0.96-1.20] | 6/3, p 0.312, 1.08 [0.96-1.20] |

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
| size2yr | squid | focal 80nAm | 45.8 | 4.2 | 0.33 | 0.08 | 0.29 |
| size2yr | squid | focal 320nAm | 16.6 | 6.3 | 0.88 | 0.12 | 0.62 |
| size2yr | squid | patch 80nAm | 49.6 | 6.9 | 0.04 | 0.04 | 0.04 |
| size2yr | squid | patch 320nAm | 28.8 | 8.5 | 0.54 | 0.08 | 0.29 |
| size2yr | opm_matched | focal 80nAm | 36.4 | 3.9 | 0.38 | 0.21 | 0.38 |
| size2yr | opm_matched | focal 320nAm | 16.1 | 5.6 | 0.88 | 0.17 | 0.67 |
| size2yr | opm_matched | patch 80nAm | 57.5 | 7.8 | 0.12 | 0.04 | 0.08 |
| size2yr | opm_matched | patch 320nAm | 22.0 | 7.9 | 0.54 | 0.17 | 0.33 |
| size2yr | opm_dense | focal 80nAm | 41.5 | 5.0 | 0.46 | 0.17 | 0.33 |
| size2yr | opm_dense | focal 320nAm | 16.8 | 6.1 | 0.92 | 0.17 | 0.62 |
| size2yr | opm_dense | patch 80nAm | 52.4 | 3.9 | 0.08 | 0.04 | 0.08 |
| size2yr | opm_dense | patch 320nAm | 18.3 | 7.9 | 0.58 | 0.12 | 0.38 |
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
| infant18mo | squid | focal 80nAm | 41.9 | 4.0 | 0.38 | 0.21 | 0.29 |
| infant18mo | squid | focal 320nAm | 14.8 | 6.0 | 0.96 | 0.29 | 0.67 |
| infant18mo | squid | patch 80nAm | 52.9 | 11.7 | 0.33 | 0.12 | 0.12 |
| infant18mo | squid | patch 320nAm | 12.7 | 8.2 | 0.88 | 0.42 | 0.58 |
| infant18mo | opm_matched | focal 80nAm | 32.8 | 4.8 | 0.46 | 0.33 | 0.29 |
| infant18mo | opm_matched | focal 320nAm | 9.5 | 4.1 | 0.96 | 0.54 | 0.71 |
| infant18mo | opm_matched | patch 80nAm | 41.1 | 9.9 | 0.29 | 0.21 | 0.17 |
| infant18mo | opm_matched | patch 320nAm | 8.4 | 8.6 | 0.92 | 0.58 | 0.54 |
| infant18mo | opm_dense | focal 80nAm | 26.4 | 3.5 | 0.46 | 0.42 | 0.42 |
| infant18mo | opm_dense | focal 320nAm | 8.0 | 4.7 | 0.96 | 0.58 | 0.62 |
| infant18mo | opm_dense | patch 80nAm | 40.0 | 6.4 | 0.33 | 0.21 | 0.21 |
| infant18mo | opm_dense | patch 320nAm | 9.8 | 8.0 | 0.96 | 0.54 | 0.71 |
| infant12mo | squid | focal 80nAm | 25.0 | 5.0 | 0.42 | 0.08 | 0.38 |
| infant12mo | squid | focal 320nAm | 13.2 | 5.5 | 0.96 | 0.17 | 0.71 |
| infant12mo | squid | patch 80nAm | 44.4 | 9.8 | 0.38 | 0.04 | 0.21 |
| infant12mo | squid | patch 320nAm | 16.7 | 6.1 | 0.96 | 0.29 | 0.67 |
| infant12mo | opm_matched | focal 80nAm | 22.3 | 6.1 | 0.46 | 0.33 | 0.42 |
| infant12mo | opm_matched | focal 320nAm | 9.5 | 5.9 | 0.96 | 0.54 | 0.67 |
| infant12mo | opm_matched | patch 80nAm | 22.7 | 6.3 | 0.38 | 0.21 | 0.29 |
| infant12mo | opm_matched | patch 320nAm | 6.5 | 7.0 | 0.96 | 0.67 | 0.67 |
| infant12mo | opm_dense | focal 80nAm | 20.3 | 5.0 | 0.50 | 0.29 | 0.42 |
| infant12mo | opm_dense | focal 320nAm | 7.4 | 5.4 | 1.00 | 0.62 | 0.62 |
| infant12mo | opm_dense | patch 80nAm | 50.2 | 7.8 | 0.38 | 0.12 | 0.25 |
| infant12mo | opm_dense | patch 320nAm | 7.0 | 6.8 | 1.00 | 0.75 | 0.67 |

Simulated IED-source recovery does not identify an epileptogenic zone or establish surgical benefit. Average templates of one database are not a population; the scaled adults are size-only controls.
