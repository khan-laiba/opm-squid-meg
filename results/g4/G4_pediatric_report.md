# G4 pediatric: IED detection and bounded localization in the fixed adult helmet (NEW)

Same framework, configuration and seeds as the adult (configs/g4_epilepsy.toml); the child anatomy, arrays and placement come from G3B (Neuromag at top contact, refitted OPM arrays). The adult rows are the frozen adult study (Neuromag at its measured head position; top contact would raise it by 5.5 mm). Thresholds are calibrated on each anatomy's own null data. p-values are uncorrected (24 paired detection comparisons per OPM array and anatomy); the location is the statistical unit.

## Strength for 50 % detection [nAm] (focal; practical detector at 1 false event/min)

| anatomy | detector | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|---|
| adult | squid/combined | 47 [36-63] | 86 [66-116] | 139 [112-204] | 302 [241-open] |
| adult | squid/grad | 55 [44-69] | 101 [74-133] | 151 [120-216] | none [290-open] |
| adult | squid/mag | 54 [45-66] | 83 [65-121] | 144 [117-207] | 311 [249-open] |
| adult | opm_matched/opm | 40 [32-58] | 84 [65-127] | 143 [113-204] | 320 [260-open] |
| adult | opm_dense/opm | 34 [28-51] | 72 [58-103] | 128 [103-182] | 254 [207-320] |
| school | squid/combined | 53 [40-76] | 85 [62-111] | 167 [132-243] | 299 [226-open] |
| school | squid/grad | 57 [47-73] | 98 [72-125] | 207 [152-295] | none [256-open] |
| school | squid/mag | 56 [47-69] | 83 [63-105] | 156 [126-254] | 291 [224-open] |
| school | opm_matched/opm | 42 [33-66] | 85 [64-106] | 175 [129-267] | 285 [243-open] |
| school | opm_dense/opm | 33 [29-44] | 85 [63-105] | 160 [124-238] | 275 [214-open] |
| size2yr | squid/combined | 61 [50-82] | 104 [78-127] | 200 [153-250] | 273 [224-open] |
| size2yr | squid/grad | 65 [53-89] | 112 [87-141] | 217 [172-287] | none [285-open] |
| size2yr | squid/mag | 59 [48-78] | 106 [80-129] | 210 [156-287] | 273 [224-open] |
| size2yr | opm_matched/opm | 38 [31-62] | 96 [71-116] | 186 [144-235] | 277 [238-open] |
| size2yr | opm_dense/opm | 37 [29-52] | 94 [68-126] | 171 [134-231] | 247 [202-open] |
| infant2yr | squid/combined | 43 [34-53] | 71 [61-95] | 143 [118-187] | 240 [202-292] |
| infant2yr | squid/grad | 51 [44-57] | 80 [64-107] | 160 [129-218] | 279 [241-open] |
| infant2yr | squid/mag | 46 [37-53] | 78 [63-103] | 139 [116-187] | 240 [203-289] |
| infant2yr | opm_matched/opm | 35 [30-45] | 70 [58-92] | 150 [117-202] | 238 [205-295] |
| infant2yr | opm_dense/opm | 30 [27-35] | 66 [58-88] | 138 [110-182] | 229 [196-279] |
| infant18mo | squid/combined | 36 [31-56] | 62 [53-76] | 124 [107-155] | none [250-open] |
| infant18mo | squid/grad | 47 [34-67] | 83 [61-118] | 144 [122-189] | none [288-open] |
| infant18mo | squid/mag | 36 [31-53] | 64 [54-95] | 126 [111-160] | none [254-open] |
| infant18mo | opm_matched/opm | 33 [26-48] | 60 [42-95] | 122 [106-151] | none [251-open] |
| infant18mo | opm_dense/opm | 27 [21-40] | 58 [43-88] | 118 [99-152] | 282 [214-open] |
| infant12mo | squid/combined | 43 [34-57] | 71 [57-98] | 131 [114-172] | 254 [220-320] |
| infant12mo | squid/grad | 52 [44-62] | 74 [62-102] | 133 [119-175] | none [256-open] |
| infant12mo | squid/mag | 41 [33-56] | 71 [60-98] | 132 [115-166] | 256 [219-open] |
| infant12mo | opm_matched/opm | 33 [30-44] | 63 [53-85] | 130 [109-170] | 256 [221-open] |
| infant12mo | opm_dense/opm | 29 [21-41] | 63 [50-93] | 119 [98-156] | 248 [215-320] |

## Held-out false events per minute at the 1-per-minute thresholds, and sensitivity at a matched held-out rate

Sensitivity for 40-nAm spikes at 10-30 mm with each detector's threshold set on the held-out null to 1 false event per minute (the frozen thresholds give 0.40-1.55 per minute over the anatomies and detectors, unequal between arrays).

| anatomy | detector | held-out rate at the frozen threshold | sensitivity at a matched 1 per minute |
|---|---|---|---|
| adult | squid/combined | 1.30 | 0.20 |
| adult | squid/grad | 1.20 | 0.15 |
| adult | squid/mag | 0.65 | 0.18 |
| adult | opm_matched/opm | 0.70 | 0.29 |
| adult | opm_dense/opm | 0.95 | 0.36 |
| school | squid/combined | 1.40 | 0.20 |
| school | squid/grad | 1.20 | 0.15 |
| school | squid/mag | 1.45 | 0.13 |
| school | opm_matched/opm | 0.50 | 0.29 |
| school | opm_dense/opm | 1.05 | 0.40 |
| size2yr | squid/combined | 0.40 | 0.14 |
| size2yr | squid/grad | 1.10 | 0.11 |
| size2yr | squid/mag | 0.55 | 0.14 |
| size2yr | opm_matched/opm | 0.65 | 0.29 |
| size2yr | opm_dense/opm | 1.20 | 0.31 |
| infant2yr | squid/combined | 0.70 | 0.26 |
| infant2yr | squid/grad | 0.60 | 0.18 |
| infant2yr | squid/mag | 1.55 | 0.22 |
| infant2yr | opm_matched/opm | 1.55 | 0.34 |
| infant2yr | opm_dense/opm | 1.30 | 0.43 |
| infant18mo | squid/combined | 0.90 | 0.37 |
| infant18mo | squid/grad | 0.75 | 0.31 |
| infant18mo | squid/mag | 1.15 | 0.36 |
| infant18mo | opm_matched/opm | 1.45 | 0.44 |
| infant18mo | opm_dense/opm | 1.10 | 0.50 |
| infant12mo | squid/combined | 1.20 | 0.28 |
| infant12mo | squid/grad | 0.95 | 0.18 |
| infant12mo | squid/mag | 0.95 | 0.30 |
| infant12mo | opm_matched/opm | 0.95 | 0.43 |
| infant12mo | opm_dense/opm | 0.85 | 0.44 |

## Paired OPM dense vs Neuromag combined (practical detector, 1 false event/min; locations favouring OPM / SQUID, sign-flip p, S50 ratio SQUID/OPM [95 % CI])

| anatomy | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|
| adult | 11/1, p 0.00391, 1.36 [1.10-1.58] | 9/1, p 0.0137, 1.20 [1.06-1.34] | 5/1, p 0.188, 1.09 [1.00-1.22] | 9/0, p 0.00391, 1.19 [1.00-open] |
| school | 13/0, p 0.000244, 1.60 [1.31-1.87] | 4/1, p 0.312, 1.01 [0.94-1.12] | 8/3, p 0.241, 1.04 [0.85-1.24] | 6/3, p 0.781, 1.08 [open-open] |
| size2yr | 17/0, p 1.53e-05, 1.65 [1.41-1.87] | 7/3, p 0.236, 1.11 [0.94-1.35] | 8/2, p 0.0488, 1.16 [1.01-1.35] | 8/1, p 0.0312, 1.11 [open-open] |
| infant2yr | 14/2, p 0.00153, 1.44 [1.20-1.65] | 6/1, p 0.219, 1.07 [0.97-1.23] | 7/3, p 0.465, 1.03 [0.93-1.15] | 7/2, p 0.148, 1.05 [0.96-1.13] |
| infant18mo | 16/0, p 3.05e-05, 1.33 [1.16-1.62] | 6/2, p 0.188, 1.05 [0.89-1.25] | 8/3, p 0.183, 1.05 [0.98-1.13] | 7/1, p 0.0469, > 1.13 (Neuromag does not reach 50 %) [open-open] |
| infant12mo | 14/1, p 0.000488, 1.46 [1.21-1.84] | 9/3, p 0.221, 1.12 [1.01-1.25] | 8/1, p 0.0312, 1.10 [1.03-1.23] | 3/2, p 0.75, 1.02 [0.93-open] |

## Localization (median error [mm]; joint detection + localization within 10 mm)

| anatomy | array | condition | dSPM error (all) | dSPM, MNE (all) | ECD error (detected) | detected | joint dSPM | joint ECD |
|---|---|---|---|---|---|---|---|---|
| adult | squid | focal 80nAm | 30.8 | 30.8 | 3.7 | 0.29 | 0.04 | 0.29 |
| adult | squid | focal 320nAm | 14.0 | 14.0 | 3.7 | 0.92 | 0.21 | 0.67 |
| adult | squid | patch 80nAm | 54.1 | 64.7 | 9.3 | 0.04 | 0.00 | 0.04 |
| adult | squid | patch 320nAm | 20.1 | 19.2 | 7.2 | 0.83 | 0.00 | 0.58 |
| adult | opm_matched | focal 80nAm | 23.2 | 21.7 | 6.2 | 0.38 | 0.17 | 0.25 |
| adult | opm_matched | focal 320nAm | 10.2 | 10.0 | 5.8 | 0.96 | 0.46 | 0.58 |
| adult | opm_matched | patch 80nAm | 46.6 | 51.6 | 27.0 | 0.12 | 0.04 | 0.04 |
| adult | opm_matched | patch 320nAm | 15.1 | 14.0 | 8.1 | 0.79 | 0.29 | 0.50 |
| adult | opm_dense | focal 80nAm | 21.3 | 15.7 | 5.1 | 0.38 | 0.12 | 0.29 |
| adult | opm_dense | focal 320nAm | 13.0 | 11.8 | 4.2 | 1.00 | 0.33 | 0.67 |
| adult | opm_dense | patch 80nAm | 63.0 | 42.7 | 22.4 | 0.08 | 0.04 | 0.04 |
| adult | opm_dense | patch 320nAm | 13.8 | 12.1 | 8.0 | 0.88 | 0.33 | 0.58 |
| adult | squid_mag | focal 80nAm | 53.3 | 50.1 | 5.3 | 0.29 | 0.04 | 0.29 |
| adult | squid_mag | focal 320nAm | 14.0 | 12.8 | 4.6 | 0.92 | 0.29 | 0.62 |
| adult | squid_mag | patch 80nAm | 66.1 | 66.5 | 9.1 | 0.04 | 0.00 | 0.04 |
| adult | squid_mag | patch 320nAm | 15.0 | 14.4 | 8.0 | 0.83 | 0.08 | 0.54 |
| adult | squid_grad | focal 80nAm | 31.3 | 26.4 | 4.7 | 0.25 | 0.04 | 0.21 |
| adult | squid_grad | focal 320nAm | 16.2 | 14.2 | 4.4 | 0.92 | 0.12 | 0.75 |
| adult | squid_grad | patch 80nAm | 93.5 | 94.5 | - | 0.00 | 0.00 | 0.00 |
| adult | squid_grad | patch 320nAm | 19.2 | 14.6 | 7.6 | 0.67 | 0.04 | 0.46 |
| school | squid | focal 80nAm | 53.4 | 51.3 | 5.5 | 0.38 | 0.00 | 0.29 |
| school | squid | focal 320nAm | 17.5 | 15.0 | 7.6 | 0.88 | 0.21 | 0.50 |
| school | squid | patch 80nAm | 55.7 | 55.7 | 6.7 | 0.08 | 0.04 | 0.08 |
| school | squid | patch 320nAm | 18.4 | 19.5 | 8.9 | 0.71 | 0.21 | 0.46 |
| school | opm_matched | focal 80nAm | 16.5 | 16.5 | 4.5 | 0.33 | 0.17 | 0.29 |
| school | opm_matched | focal 320nAm | 10.5 | 11.2 | 7.6 | 0.88 | 0.38 | 0.58 |
| school | opm_matched | patch 80nAm | 49.9 | 51.0 | 8.8 | 0.08 | 0.04 | 0.04 |
| school | opm_matched | patch 320nAm | 10.8 | 11.2 | 8.1 | 0.75 | 0.33 | 0.54 |
| school | opm_dense | focal 80nAm | 21.5 | 17.4 | 4.2 | 0.42 | 0.08 | 0.33 |
| school | opm_dense | focal 320nAm | 10.1 | 10.1 | 6.4 | 0.96 | 0.50 | 0.62 |
| school | opm_dense | patch 80nAm | 55.2 | 51.4 | 10.2 | 0.17 | 0.04 | 0.08 |
| school | opm_dense | patch 320nAm | 12.8 | 11.1 | 8.2 | 0.79 | 0.29 | 0.54 |
| school | squid_mag | focal 80nAm | 20.9 | 22.1 | 5.2 | 0.38 | 0.17 | 0.29 |
| school | squid_mag | focal 320nAm | 13.8 | 11.1 | 8.5 | 0.83 | 0.25 | 0.46 |
| school | squid_mag | patch 80nAm | 51.3 | 55.4 | 6.0 | 0.08 | 0.04 | 0.08 |
| school | squid_mag | patch 320nAm | 11.7 | 11.2 | 9.7 | 0.75 | 0.38 | 0.46 |
| school | squid_grad | focal 80nAm | 50.1 | 55.6 | 5.5 | 0.33 | 0.08 | 0.25 |
| school | squid_grad | focal 320nAm | 15.6 | 13.6 | 8.0 | 0.83 | 0.33 | 0.54 |
| school | squid_grad | patch 80nAm | 57.8 | 55.2 | 6.9 | 0.08 | 0.04 | 0.08 |
| school | squid_grad | patch 320nAm | 16.8 | 19.5 | 9.1 | 0.67 | 0.25 | 0.38 |
| size2yr | squid | focal 80nAm | 38.5 | 35.9 | 4.2 | 0.33 | 0.08 | 0.29 |
| size2yr | squid | focal 320nAm | 16.6 | 16.6 | 6.0 | 0.88 | 0.12 | 0.62 |
| size2yr | squid | patch 80nAm | 53.1 | 53.1 | 6.8 | 0.04 | 0.04 | 0.04 |
| size2yr | squid | patch 320nAm | 21.0 | 20.7 | 7.5 | 0.50 | 0.08 | 0.29 |
| size2yr | opm_matched | focal 80nAm | 36.4 | 36.0 | 3.9 | 0.38 | 0.21 | 0.38 |
| size2yr | opm_matched | focal 320nAm | 16.1 | 12.4 | 5.6 | 0.88 | 0.17 | 0.62 |
| size2yr | opm_matched | patch 80nAm | 57.5 | 54.2 | 8.2 | 0.12 | 0.04 | 0.08 |
| size2yr | opm_matched | patch 320nAm | 22.0 | 17.6 | 7.9 | 0.54 | 0.17 | 0.33 |
| size2yr | opm_dense | focal 80nAm | 38.1 | 40.8 | 4.5 | 0.46 | 0.25 | 0.38 |
| size2yr | opm_dense | focal 320nAm | 16.1 | 12.1 | 6.1 | 0.92 | 0.21 | 0.62 |
| size2yr | opm_dense | patch 80nAm | 47.3 | 51.0 | 6.2 | 0.12 | 0.04 | 0.08 |
| size2yr | opm_dense | patch 320nAm | 14.3 | 18.8 | 7.3 | 0.54 | 0.12 | 0.38 |
| size2yr | squid_mag | focal 80nAm | 38.9 | 40.8 | 4.4 | 0.33 | 0.08 | 0.29 |
| size2yr | squid_mag | focal 320nAm | 14.3 | 13.2 | 6.2 | 0.88 | 0.21 | 0.54 |
| size2yr | squid_mag | patch 80nAm | 48.1 | 44.8 | 11.8 | 0.08 | 0.04 | 0.04 |
| size2yr | squid_mag | patch 320nAm | 25.0 | 25.0 | 7.4 | 0.50 | 0.08 | 0.29 |
| size2yr | squid_grad | focal 80nAm | 31.5 | 38.1 | 7.6 | 0.33 | 0.04 | 0.17 |
| size2yr | squid_grad | focal 320nAm | 15.9 | 14.6 | 6.8 | 0.79 | 0.17 | 0.50 |
| size2yr | squid_grad | patch 80nAm | 46.3 | 48.8 | 6.5 | 0.04 | 0.00 | 0.04 |
| size2yr | squid_grad | patch 320nAm | 26.7 | 34.0 | 6.0 | 0.42 | 0.04 | 0.38 |
| infant2yr | squid | focal 80nAm | 20.1 | 34.2 | 6.8 | 0.42 | 0.12 | 0.38 |
| infant2yr | squid | focal 320nAm | 14.8 | 13.6 | 4.9 | 0.83 | 0.29 | 0.67 |
| infant2yr | squid | patch 80nAm | 41.3 | 41.3 | 7.1 | 0.33 | 0.04 | 0.25 |
| infant2yr | squid | patch 320nAm | 13.6 | 11.4 | 6.2 | 0.88 | 0.42 | 0.75 |
| infant2yr | opm_matched | focal 80nAm | 22.3 | 24.4 | 6.5 | 0.50 | 0.33 | 0.38 |
| infant2yr | opm_matched | focal 320nAm | 8.9 | 7.6 | 4.7 | 0.92 | 0.54 | 0.75 |
| infant2yr | opm_matched | patch 80nAm | 37.0 | 42.6 | 4.6 | 0.29 | 0.12 | 0.25 |
| infant2yr | opm_matched | patch 320nAm | 10.0 | 9.8 | 4.9 | 0.96 | 0.50 | 0.67 |
| infant2yr | opm_dense | focal 80nAm | 16.6 | 14.8 | 5.2 | 0.50 | 0.29 | 0.38 |
| infant2yr | opm_dense | focal 320nAm | 7.5 | 6.4 | 5.2 | 0.92 | 0.54 | 0.58 |
| infant2yr | opm_dense | patch 80nAm | 39.0 | 40.5 | 6.2 | 0.33 | 0.12 | 0.25 |
| infant2yr | opm_dense | patch 320nAm | 10.1 | 8.9 | 7.2 | 1.00 | 0.50 | 0.71 |
| infant2yr | squid_mag | focal 80nAm | 18.3 | 24.4 | 6.5 | 0.38 | 0.12 | 0.29 |
| infant2yr | squid_mag | focal 320nAm | 13.3 | 12.2 | 4.1 | 0.79 | 0.38 | 0.67 |
| infant2yr | squid_mag | patch 80nAm | 37.5 | 35.9 | 7.9 | 0.33 | 0.00 | 0.21 |
| infant2yr | squid_mag | patch 320nAm | 11.0 | 10.0 | 5.5 | 0.88 | 0.42 | 0.79 |
| infant2yr | squid_grad | focal 80nAm | 29.2 | 20.2 | 6.1 | 0.38 | 0.17 | 0.25 |
| infant2yr | squid_grad | focal 320nAm | 15.9 | 15.9 | 4.5 | 0.79 | 0.29 | 0.67 |
| infant2yr | squid_grad | patch 80nAm | 57.7 | 57.8 | 8.8 | 0.17 | 0.04 | 0.12 |
| infant2yr | squid_grad | patch 320nAm | 13.7 | 13.7 | 6.6 | 0.88 | 0.33 | 0.67 |
| infant18mo | squid | focal 80nAm | 38.0 | 45.7 | 4.0 | 0.38 | 0.21 | 0.29 |
| infant18mo | squid | focal 320nAm | 14.8 | 13.2 | 5.9 | 0.96 | 0.29 | 0.71 |
| infant18mo | squid | patch 80nAm | 56.4 | 56.4 | 12.4 | 0.33 | 0.12 | 0.12 |
| infant18mo | squid | patch 320nAm | 12.7 | 12.7 | 8.2 | 0.88 | 0.42 | 0.54 |
| infant18mo | opm_matched | focal 80nAm | 32.8 | 34.9 | 4.8 | 0.46 | 0.33 | 0.29 |
| infant18mo | opm_matched | focal 320nAm | 9.5 | 8.3 | 4.1 | 0.96 | 0.54 | 0.67 |
| infant18mo | opm_matched | patch 80nAm | 41.1 | 30.5 | 9.8 | 0.29 | 0.21 | 0.17 |
| infant18mo | opm_matched | patch 320nAm | 8.4 | 8.2 | 9.2 | 0.92 | 0.58 | 0.54 |
| infant18mo | opm_dense | focal 80nAm | 26.4 | 35.5 | 4.7 | 0.46 | 0.42 | 0.42 |
| infant18mo | opm_dense | focal 320nAm | 8.0 | 8.0 | 4.8 | 0.96 | 0.58 | 0.62 |
| infant18mo | opm_dense | patch 80nAm | 40.0 | 36.8 | 7.7 | 0.33 | 0.21 | 0.17 |
| infant18mo | opm_dense | patch 320nAm | 9.2 | 10.0 | 7.4 | 0.96 | 0.54 | 0.71 |
| infant18mo | squid_mag | focal 80nAm | 22.0 | 26.7 | 5.5 | 0.42 | 0.17 | 0.29 |
| infant18mo | squid_mag | focal 320nAm | 13.0 | 13.0 | 5.5 | 0.92 | 0.33 | 0.71 |
| infant18mo | squid_mag | patch 80nAm | 53.3 | 49.7 | 12.3 | 0.33 | 0.08 | 0.12 |
| infant18mo | squid_mag | patch 320nAm | 10.9 | 13.7 | 7.8 | 0.88 | 0.46 | 0.58 |
| infant18mo | squid_grad | focal 80nAm | 48.7 | 48.2 | 4.5 | 0.38 | 0.17 | 0.29 |
| infant18mo | squid_grad | focal 320nAm | 14.3 | 13.4 | 5.9 | 0.83 | 0.25 | 0.67 |
| infant18mo | squid_grad | patch 80nAm | 58.1 | 61.6 | 9.5 | 0.33 | 0.08 | 0.17 |
| infant18mo | squid_grad | patch 320nAm | 17.0 | 16.1 | 7.5 | 0.79 | 0.25 | 0.62 |
| infant12mo | squid | focal 80nAm | 39.6 | 43.1 | 7.9 | 0.50 | 0.08 | 0.38 |
| infant12mo | squid | focal 320nAm | 12.0 | 12.4 | 5.6 | 0.96 | 0.42 | 0.83 |
| infant12mo | squid | patch 80nAm | 59.7 | 59.4 | 7.8 | 0.33 | 0.00 | 0.21 |
| infant12mo | squid | patch 320nAm | 12.8 | 12.5 | 6.3 | 0.96 | 0.38 | 0.67 |
| infant12mo | opm_matched | focal 80nAm | 11.9 | 10.3 | 6.2 | 0.54 | 0.38 | 0.50 |
| infant12mo | opm_matched | focal 320nAm | 10.0 | 7.4 | 5.3 | 0.96 | 0.50 | 0.75 |
| infant12mo | opm_matched | patch 80nAm | 16.7 | 20.8 | 4.7 | 0.38 | 0.29 | 0.29 |
| infant12mo | opm_matched | patch 320nAm | 8.2 | 10.0 | 5.2 | 0.96 | 0.58 | 0.75 |
| infant12mo | opm_dense | focal 80nAm | 12.6 | 11.2 | 5.5 | 0.54 | 0.42 | 0.38 |
| infant12mo | opm_dense | focal 320nAm | 8.7 | 7.4 | 5.2 | 1.00 | 0.58 | 0.79 |
| infant12mo | opm_dense | patch 80nAm | 24.8 | 23.0 | 7.7 | 0.46 | 0.25 | 0.25 |
| infant12mo | opm_dense | patch 320nAm | 7.4 | 8.2 | 7.7 | 0.96 | 0.54 | 0.67 |
| infant12mo | squid_mag | focal 80nAm | 24.6 | 26.0 | 8.0 | 0.50 | 0.12 | 0.33 |
| infant12mo | squid_mag | focal 320nAm | 11.9 | 11.7 | 5.2 | 0.96 | 0.33 | 0.75 |
| infant12mo | squid_mag | patch 80nAm | 19.5 | 21.2 | 7.7 | 0.33 | 0.12 | 0.21 |
| infant12mo | squid_mag | patch 320nAm | 11.5 | 9.9 | 7.2 | 0.96 | 0.42 | 0.71 |
| infant12mo | squid_grad | focal 80nAm | 35.5 | 41.8 | 10.6 | 0.46 | 0.08 | 0.21 |
| infant12mo | squid_grad | focal 320nAm | 15.8 | 12.1 | 5.3 | 0.96 | 0.33 | 0.79 |
| infant12mo | squid_grad | patch 80nAm | 42.3 | 50.8 | 7.2 | 0.25 | 0.04 | 0.17 |
| infant12mo | squid_grad | patch 320nAm | 16.4 | 15.5 | 7.5 | 0.96 | 0.25 | 0.58 |

Localization detectors (1 false event per minute on 10 min of null data) on independent held-out null data: false events per minute [exact 95 % interval]:

| anatomy | squid | opm_matched | opm_dense | squid_mag | squid_grad |
|---|---|---|---|---|---|
| adult | 1.30 [0.69-2.22] | 1.60 [0.91-2.60] | 1.10 [0.55-1.97] | 1.30 [0.69-2.22] | 1.10 [0.55-1.97] |
| school | 0.70 [0.28-1.44] | 0.30 [0.06-0.88] | 1.20 [0.62-2.10] | 0.50 [0.16-1.17] | 1.00 [0.48-1.84] |
| size2yr | 1.00 [0.48-1.84] | 1.20 [0.62-2.10] | 1.00 [0.48-1.84] | 1.10 [0.55-1.97] | 1.20 [0.62-2.10] |
| infant2yr | 1.00 [0.48-1.84] | 1.30 [0.69-2.22] | 1.60 [0.91-2.60] | 0.70 [0.28-1.44] | 1.50 [0.84-2.47] |
| infant18mo | 1.30 [0.69-2.22] | 1.50 [0.84-2.47] | 0.80 [0.35-1.58] | 1.10 [0.55-1.97] | 1.30 [0.69-2.22] |
| infant12mo | 1.50 [0.84-2.47] | 0.70 [0.28-1.44] | 1.40 [0.77-2.35] | 1.80 [1.07-2.84] | 1.10 [0.55-1.97] |

Simulated IED-source recovery does not identify an epileptogenic zone or establish surgical benefit. Average templates of one database are not a population; the scaled adults are size-only controls.
