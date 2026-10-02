# G4 pediatric: IED detection and bounded localization in the fixed adult helmet (NEW)

Same framework, configuration and seeds as the adult (configs/g4_epilepsy.toml); the child anatomy, arrays and placement come from G3B (Neuromag at top contact, refitted OPM arrays). The adult rows are the frozen adult study (Neuromag at its measured head position; top contact would raise it by 5.5 mm). Thresholds are calibrated on each anatomy's own null data. p-values are uncorrected (24 paired detection comparisons per OPM array and anatomy); the location is the statistical unit.

## Strength for 50 % detection [nAm] (focal; practical detector at 1 false event/min)

| anatomy | detector | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|---|
| adult | squid/combined | 48 [37-63] | 85 [68-110] | 136 [110-194] | 263 [open] |
| adult | squid/grad | 57 [46-71] | 101 [75-133] | 143 [117-211] | none [open] |
| adult | squid/mag | 52 [40-66] | 86 [66-113] | 136 [111-189] | 275 [open] |
| adult | opm_matched/opm | 44 [33-61] | 87 [64-123] | 146 [112-210] | 320 [open] |
| adult | opm_dense/opm | 37 [30-55] | 80 [61-109] | 135 [110-185] | 268 [open] |
| school | squid/combined | 50 [39-67] | 92 [69-113] | 167 [129-244] | 271 [open] |
| school | squid/grad | 58 [50-71] | 94 [72-116] | 202 [144-300] | none [open] |
| school | squid/mag | 54 [44-72] | 92 [70-108] | 179 [131-247] | 271 [open] |
| school | opm_matched/opm | 44 [33-71] | 89 [69-106] | 185 [136-258] | 311 [open] |
| school | opm_dense/opm | 36 [30-50] | 86 [69-104] | 151 [121-231] | 284 [open] |
| size2yr | squid/combined | 59 [49-77] | 106 [78-135] | 200 [153-250] | 273 [open] |
| size2yr | squid/grad | 65 [53-89] | 110 [87-135] | 217 [172-287] | none [open] |
| size2yr | squid/mag | 59 [48-78] | 107 [80-134] | 204 [152-269] | 273 [open] |
| size2yr | opm_matched/opm | 38 [31-62] | 96 [71-116] | 187 [144-239] | 279 [open] |
| size2yr | opm_dense/opm | 37 [29-51] | 89 [65-120] | 177 [134-232] | 243 [open] |
| infant2yr | squid/combined | 43 [34-53] | 71 [61-95] | 143 [118-187] | 240 [202-292] |
| infant2yr | squid/grad | 51 [44-57] | 80 [64-107] | 160 [129-218] | 279 [open] |
| infant2yr | squid/mag | 46 [37-53] | 78 [63-103] | 139 [116-187] | 240 [203-289] |
| infant2yr | opm_matched/opm | 35 [30-45] | 70 [58-92] | 150 [117-202] | 238 [205-295] |
| infant2yr | opm_dense/opm | 30 [27-35] | 66 [58-88] | 138 [110-182] | 229 [196-279] |
| infant18mo | squid/combined | 36 [31-56] | 62 [53-76] | 124 [107-155] | none [open] |
| infant18mo | squid/grad | 47 [34-67] | 83 [61-118] | 144 [122-189] | none [open] |
| infant18mo | squid/mag | 36 [31-53] | 64 [54-95] | 126 [111-160] | none [open] |
| infant18mo | opm_matched/opm | 33 [26-48] | 60 [42-95] | 122 [106-151] | none [open] |
| infant18mo | opm_dense/opm | 27 [21-40] | 58 [43-88] | 118 [99-152] | 282 [open] |
| infant12mo | squid/combined | 43 [34-57] | 71 [57-98] | 131 [114-172] | 254 [220-320] |
| infant12mo | squid/grad | 52 [44-62] | 74 [62-102] | 133 [119-175] | none [open] |
| infant12mo | squid/mag | 41 [33-56] | 71 [60-98] | 132 [115-166] | 256 [open] |
| infant12mo | opm_matched/opm | 33 [30-44] | 63 [53-85] | 130 [109-170] | 256 [open] |
| infant12mo | opm_dense/opm | 29 [21-41] | 63 [50-93] | 119 [98-156] | 248 [215-320] |

## Held-out false events per minute at the 1-per-minute thresholds, and sensitivity at a matched held-out rate

Sensitivity for 40-nAm spikes at 10-30 mm with each detector's threshold set on the held-out null to 1 false event per minute (the frozen thresholds give 0.5-1.7 per minute, unequal between arrays).

| anatomy | detector | held-out rate at the frozen threshold | sensitivity at a matched 1 per minute |
|---|---|---|---|
| adult | squid/combined | 1.00 | 0.21 |
| adult | squid/grad | 0.80 | 0.12 |
| adult | squid/mag | 1.15 | 0.17 |
| adult | opm_matched/opm | 0.65 | 0.28 |
| adult | opm_dense/opm | 0.80 | 0.33 |
| school | squid/combined | 1.00 | 0.20 |
| school | squid/grad | 0.95 | 0.12 |
| school | squid/mag | 1.20 | 0.17 |
| school | opm_matched/opm | 0.70 | 0.24 |
| school | opm_dense/opm | 0.40 | 0.31 |
| size2yr | squid/combined | 0.45 | 0.13 |
| size2yr | squid/grad | 1.10 | 0.11 |
| size2yr | squid/mag | 0.60 | 0.13 |
| size2yr | opm_matched/opm | 0.70 | 0.29 |
| size2yr | opm_dense/opm | 1.25 | 0.31 |
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
| infant12mo | squid/combined | 1.25 | 0.28 |
| infant12mo | squid/grad | 0.95 | 0.18 |
| infant12mo | squid/mag | 0.95 | 0.30 |
| infant12mo | opm_matched/opm | 0.95 | 0.43 |
| infant12mo | opm_dense/opm | 0.85 | 0.44 |

## Paired OPM dense vs Neuromag combined (practical detector, 1 false event/min; locations favouring OPM / SQUID, sign-flip p, S50 ratio SQUID/OPM [95 % CI])

| anatomy | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|
| adult | 8/2, p 0.0371, 1.29 [1.03-1.51] | 5/1, p 0.188, 1.07 [0.95-1.17] | 4/0, p 0.125, 1.01 [0.93-1.08] | 0/2, p 0.5, 0.98 [open-open] |
| school | 11/1, p 0.00293, 1.40 [1.16-1.60] | 6/2, p 0.344, 1.07 [0.97-1.17] | 6/3, p 0.398, 1.10 [0.96-1.26] | 1/4, p 0.375, 0.96 [open-open] |
| size2yr | 16/0, p 3.05e-05, 1.61 [1.33-1.85] | 9/2, p 0.0449, 1.19 [0.98-1.46] | 9/2, p 0.0449, 1.13 [0.99-1.35] | 9/1, p 0.0195, 1.12 [1.00-open] |
| infant2yr | 14/2, p 0.00153, 1.44 [1.20-1.65] | 6/1, p 0.219, 1.07 [0.97-1.23] | 7/3, p 0.465, 1.03 [0.93-1.15] | 7/2, p 0.148, 1.05 [0.96-1.13] |
| infant18mo | 16/0, p 3.05e-05, 1.33 [1.16-1.62] | 6/2, p 0.188, 1.05 [0.89-1.25] | 8/3, p 0.183, 1.05 [0.98-1.13] | 7/1, p 0.0469, > 1.13 (Neuromag does not reach 50 %) [open-open] |
| infant12mo | 14/1, p 0.000488, 1.46 [1.21-1.84] | 9/3, p 0.221, 1.12 [1.01-1.25] | 8/1, p 0.0312, 1.10 [1.03-1.23] | 3/2, p 0.75, 1.02 [0.93-open] |

## Localization (median error [mm]; joint detection + localization within 10 mm)

| anatomy | array | condition | dSPM error (all) | ECD error (detected) | detected | joint dSPM | joint ECD |
|---|---|---|---|---|---|---|---|
| adult | squid | focal 80nAm | 41.1 | 4.4 | 0.33 | 0.04 | 0.29 |
| adult | squid | focal 320nAm | 13.0 | 4.7 | 0.96 | 0.25 | 0.75 |
| adult | squid | patch 80nAm | 73.8 | 16.3 | 0.08 | 0.00 | 0.04 |
| adult | squid | patch 320nAm | 14.5 | 7.0 | 0.67 | 0.21 | 0.54 |
| adult | opm_matched | focal 80nAm | 29.1 | 4.8 | 0.29 | 0.12 | 0.21 |
| adult | opm_matched | focal 320nAm | 11.1 | 4.2 | 0.79 | 0.38 | 0.67 |
| adult | opm_matched | patch 80nAm | 55.2 | 8.4 | 0.08 | 0.00 | 0.04 |
| adult | opm_matched | patch 320nAm | 15.3 | 8.2 | 0.71 | 0.29 | 0.50 |
| adult | opm_dense | focal 80nAm | 17.3 | 5.3 | 0.38 | 0.17 | 0.33 |
| adult | opm_dense | focal 320nAm | 11.6 | 4.6 | 0.96 | 0.42 | 0.71 |
| adult | opm_dense | patch 80nAm | 71.1 | 20.6 | 0.12 | 0.00 | 0.04 |
| adult | opm_dense | patch 320nAm | 13.4 | 6.2 | 0.75 | 0.38 | 0.62 |
| school | squid | focal 80nAm | 33.9 | 9.2 | 0.46 | 0.00 | 0.25 |
| school | squid | focal 320nAm | 14.0 | 7.9 | 0.96 | 0.29 | 0.50 |
| school | squid | patch 80nAm | 63.3 | 57.4 | 0.04 | 0.00 | 0.00 |
| school | squid | patch 320nAm | 20.5 | 7.4 | 0.71 | 0.17 | 0.50 |
| school | opm_matched | focal 80nAm | 33.1 | 5.2 | 0.42 | 0.25 | 0.33 |
| school | opm_matched | focal 320nAm | 11.3 | 7.6 | 0.96 | 0.33 | 0.62 |
| school | opm_matched | patch 80nAm | 54.9 | - | 0.00 | 0.00 | 0.00 |
| school | opm_matched | patch 320nAm | 10.5 | 7.6 | 0.71 | 0.38 | 0.54 |
| school | opm_dense | focal 80nAm | 27.9 | 5.4 | 0.46 | 0.08 | 0.29 |
| school | opm_dense | focal 320nAm | 12.6 | 6.0 | 0.96 | 0.33 | 0.58 |
| school | opm_dense | patch 80nAm | 58.7 | - | 0.00 | 0.00 | 0.00 |
| school | opm_dense | patch 320nAm | 13.3 | 7.4 | 0.71 | 0.38 | 0.54 |
| size2yr | squid | focal 80nAm | 38.5 | 4.2 | 0.33 | 0.08 | 0.29 |
| size2yr | squid | focal 320nAm | 16.6 | 6.0 | 0.88 | 0.12 | 0.62 |
| size2yr | squid | patch 80nAm | 53.1 | 6.7 | 0.04 | 0.04 | 0.04 |
| size2yr | squid | patch 320nAm | 21.0 | 7.3 | 0.50 | 0.08 | 0.29 |
| size2yr | opm_matched | focal 80nAm | 36.4 | 3.9 | 0.38 | 0.21 | 0.38 |
| size2yr | opm_matched | focal 320nAm | 16.1 | 5.5 | 0.88 | 0.17 | 0.62 |
| size2yr | opm_matched | patch 80nAm | 57.5 | 8.0 | 0.12 | 0.04 | 0.08 |
| size2yr | opm_matched | patch 320nAm | 22.0 | 7.9 | 0.54 | 0.17 | 0.33 |
| size2yr | opm_dense | focal 80nAm | 44.1 | 6.2 | 0.42 | 0.17 | 0.25 |
| size2yr | opm_dense | focal 320nAm | 12.9 | 5.5 | 0.92 | 0.29 | 0.62 |
| size2yr | opm_dense | patch 80nAm | 47.8 | 6.6 | 0.12 | 0.04 | 0.08 |
| size2yr | opm_dense | patch 320nAm | 14.3 | 6.5 | 0.50 | 0.17 | 0.33 |
| infant2yr | squid | focal 80nAm | 20.1 | 6.8 | 0.42 | 0.12 | 0.33 |
| infant2yr | squid | focal 320nAm | 14.8 | 4.9 | 0.83 | 0.29 | 0.67 |
| infant2yr | squid | patch 80nAm | 41.3 | 7.6 | 0.33 | 0.04 | 0.21 |
| infant2yr | squid | patch 320nAm | 13.6 | 5.9 | 0.88 | 0.42 | 0.75 |
| infant2yr | opm_matched | focal 80nAm | 22.3 | 6.5 | 0.50 | 0.33 | 0.33 |
| infant2yr | opm_matched | focal 320nAm | 8.9 | 4.8 | 0.92 | 0.54 | 0.75 |
| infant2yr | opm_matched | patch 80nAm | 37.0 | 4.6 | 0.29 | 0.12 | 0.25 |
| infant2yr | opm_matched | patch 320nAm | 10.0 | 5.0 | 0.96 | 0.50 | 0.58 |
| infant2yr | opm_dense | focal 80nAm | 16.6 | 5.2 | 0.50 | 0.29 | 0.38 |
| infant2yr | opm_dense | focal 320nAm | 7.5 | 5.2 | 0.92 | 0.54 | 0.58 |
| infant2yr | opm_dense | patch 80nAm | 39.0 | 6.2 | 0.33 | 0.12 | 0.25 |
| infant2yr | opm_dense | patch 320nAm | 10.1 | 7.2 | 1.00 | 0.50 | 0.67 |
| infant18mo | squid | focal 80nAm | 38.0 | 4.0 | 0.38 | 0.21 | 0.29 |
| infant18mo | squid | focal 320nAm | 14.8 | 6.0 | 0.96 | 0.29 | 0.67 |
| infant18mo | squid | patch 80nAm | 56.4 | 12.4 | 0.33 | 0.12 | 0.12 |
| infant18mo | squid | patch 320nAm | 12.7 | 8.2 | 0.88 | 0.42 | 0.58 |
| infant18mo | opm_matched | focal 80nAm | 32.8 | 4.8 | 0.46 | 0.33 | 0.29 |
| infant18mo | opm_matched | focal 320nAm | 9.5 | 4.1 | 0.96 | 0.54 | 0.67 |
| infant18mo | opm_matched | patch 80nAm | 41.1 | 9.9 | 0.29 | 0.21 | 0.17 |
| infant18mo | opm_matched | patch 320nAm | 8.4 | 9.1 | 0.92 | 0.58 | 0.54 |
| infant18mo | opm_dense | focal 80nAm | 26.4 | 3.6 | 0.46 | 0.42 | 0.42 |
| infant18mo | opm_dense | focal 320nAm | 8.0 | 4.8 | 0.96 | 0.58 | 0.62 |
| infant18mo | opm_dense | patch 80nAm | 40.0 | 6.4 | 0.33 | 0.21 | 0.21 |
| infant18mo | opm_dense | patch 320nAm | 9.2 | 7.5 | 0.96 | 0.54 | 0.71 |
| infant12mo | squid | focal 80nAm | 39.6 | 8.0 | 0.50 | 0.08 | 0.38 |
| infant12mo | squid | focal 320nAm | 12.0 | 5.6 | 0.96 | 0.42 | 0.79 |
| infant12mo | squid | patch 80nAm | 59.7 | 7.8 | 0.33 | 0.00 | 0.21 |
| infant12mo | squid | patch 320nAm | 12.8 | 6.3 | 0.96 | 0.38 | 0.67 |
| infant12mo | opm_matched | focal 80nAm | 11.9 | 6.2 | 0.54 | 0.38 | 0.50 |
| infant12mo | opm_matched | focal 320nAm | 10.0 | 5.4 | 0.96 | 0.50 | 0.75 |
| infant12mo | opm_matched | patch 80nAm | 16.7 | 4.7 | 0.38 | 0.29 | 0.29 |
| infant12mo | opm_matched | patch 320nAm | 8.2 | 6.2 | 0.96 | 0.58 | 0.75 |
| infant12mo | opm_dense | focal 80nAm | 12.6 | 5.5 | 0.54 | 0.42 | 0.38 |
| infant12mo | opm_dense | focal 320nAm | 8.7 | 5.2 | 1.00 | 0.58 | 0.79 |
| infant12mo | opm_dense | patch 80nAm | 24.8 | 7.7 | 0.46 | 0.25 | 0.25 |
| infant12mo | opm_dense | patch 320nAm | 7.4 | 7.5 | 0.96 | 0.54 | 0.71 |

Simulated IED-source recovery does not identify an epileptogenic zone or establish surgical benefit. Average templates of one database are not a population; the scaled adults are size-only controls.
