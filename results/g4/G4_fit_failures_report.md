# G4 localization: failed fits under declared criteria

Every event has a dipole (`mne.fit_dipole` never failed). Counted here, from the stored per-event tables: dipole (ECD) or dSPM-peak errors above 30 mm (gross), among all events and among the detected ones, and dipole confidence volumes above 10 cm^3 (unconstrained). The limits are operational choices. Conditions: 24 events each (focal and 10-mm patches at 80 and 320 nAm).

| anatomy | array | condition | detected | ECD > 30 mm (all / detected) | dSPM > 30 mm (all / detected) | ECD volume > 10 cm^3 (all / detected) |
|---|---|---|---|---|---|---|
| adult | opm_dense | focal 80nAm | 9/24 | 10 / 0 | 10 / 0 | 10 / 1 |
| adult | opm_dense | focal 320nAm | 23/24 | 3 / 2 | 5 / 4 | 2 / 1 |
| adult | opm_dense | patch 80nAm | 3/24 | 15 / 0 | 20 / 0 | 13 / 1 |
| adult | opm_dense | patch 320nAm | 18/24 | 4 / 0 | 6 / 2 | 3 / 0 |
| adult | opm_matched | focal 80nAm | 7/24 | 12 / 0 | 12 / 0 | 19 / 2 |
| adult | opm_matched | focal 320nAm | 19/24 | 4 / 1 | 4 / 1 | 6 / 2 |
| adult | opm_matched | patch 80nAm | 2/24 | 18 / 0 | 17 / 0 | 22 / 1 |
| adult | opm_matched | patch 320nAm | 17/24 | 5 / 0 | 7 / 1 | 6 / 1 |
| adult | squid | focal 80nAm | 8/24 | 12 / 0 | 12 / 1 | 16 / 1 |
| adult | squid | focal 320nAm | 23/24 | 4 / 3 | 4 / 3 | 5 / 4 |
| adult | squid | patch 80nAm | 2/24 | 18 / 0 | 20 / 1 | 20 / 1 |
| adult | squid | patch 320nAm | 16/24 | 5 / 0 | 7 / 2 | 6 / 0 |
| school | opm_dense | focal 80nAm | 11/24 | 7 / 0 | 12 / 1 | 10 / 2 |
| school | opm_dense | focal 320nAm | 23/24 | 5 / 4 | 4 / 3 | 6 / 5 |
| school | opm_dense | patch 80nAm | 0/24 | 17 / 0 | 16 / 0 | 12 / 0 |
| school | opm_dense | patch 320nAm | 17/24 | 2 / 0 | 6 / 2 | 5 / 1 |
| school | opm_matched | focal 80nAm | 10/24 | 11 / 0 | 12 / 1 | 13 / 2 |
| school | opm_matched | focal 320nAm | 23/24 | 4 / 3 | 5 / 4 | 6 / 5 |
| school | opm_matched | patch 80nAm | 0/24 | 18 / 0 | 17 / 0 | 22 / 0 |
| school | opm_matched | patch 320nAm | 17/24 | 4 / 0 | 3 / 0 | 8 / 2 |
| school | squid | focal 80nAm | 11/24 | 9 / 2 | 12 / 4 | 15 / 2 |
| school | squid | focal 320nAm | 23/24 | 4 / 3 | 6 / 5 | 6 / 5 |
| school | squid | patch 80nAm | 1/24 | 17 / 1 | 22 / 1 | 21 / 1 |
| school | squid | patch 320nAm | 17/24 | 6 / 1 | 11 / 4 | 8 / 2 |
| size2yr | opm_dense | focal 80nAm | 10/24 | 14 / 1 | 14 / 1 | 13 / 1 |
| size2yr | opm_dense | focal 320nAm | 22/24 | 1 / 1 | 2 / 2 | 4 / 2 |
| size2yr | opm_dense | patch 80nAm | 3/24 | 19 / 1 | 19 / 0 | 15 / 0 |
| size2yr | opm_dense | patch 320nAm | 12/24 | 10 / 0 | 8 / 0 | 7 / 1 |
| size2yr | opm_matched | focal 80nAm | 9/24 | 11 / 0 | 13 / 0 | 14 / 1 |
| size2yr | opm_matched | focal 320nAm | 21/24 | 2 / 2 | 5 / 3 | 4 / 2 |
| size2yr | opm_matched | patch 80nAm | 3/24 | 18 / 1 | 20 / 1 | 17 / 1 |
| size2yr | opm_matched | patch 320nAm | 13/24 | 9 / 2 | 10 / 1 | 11 / 1 |
| size2yr | squid | focal 80nAm | 8/24 | 17 / 1 | 13 / 0 | 17 / 1 |
| size2yr | squid | focal 320nAm | 21/24 | 2 / 2 | 8 / 6 | 7 / 4 |
| size2yr | squid | patch 80nAm | 1/24 | 19 / 0 | 17 / 0 | 21 / 0 |
| size2yr | squid | patch 320nAm | 12/24 | 9 / 0 | 11 / 1 | 14 / 2 |
| infant2yr | opm_dense | focal 80nAm | 12/24 | 9 / 0 | 9 / 0 | 6 / 0 |
| infant2yr | opm_dense | focal 320nAm | 22/24 | 5 / 4 | 3 / 3 | 2 / 1 |
| infant2yr | opm_dense | patch 80nAm | 8/24 | 12 / 0 | 12 / 1 | 9 / 1 |
| infant2yr | opm_dense | patch 320nAm | 24/24 | 3 / 3 | 2 / 2 | 2 / 2 |
| infant2yr | opm_matched | focal 80nAm | 12/24 | 11 / 0 | 11 / 1 | 12 / 2 |
| infant2yr | opm_matched | focal 320nAm | 22/24 | 5 / 4 | 4 / 3 | 3 / 2 |
| infant2yr | opm_matched | patch 80nAm | 7/24 | 12 / 0 | 13 / 0 | 15 / 2 |
| infant2yr | opm_matched | patch 320nAm | 23/24 | 3 / 2 | 3 / 2 | 4 / 3 |
| infant2yr | squid | focal 80nAm | 10/24 | 10 / 0 | 11 / 0 | 14 / 2 |
| infant2yr | squid | focal 320nAm | 20/24 | 5 / 2 | 5 / 2 | 5 / 2 |
| infant2yr | squid | patch 80nAm | 8/24 | 14 / 1 | 16 / 3 | 15 / 2 |
| infant2yr | squid | patch 320nAm | 21/24 | 4 / 1 | 4 / 2 | 5 / 2 |
| infant18mo | opm_dense | focal 80nAm | 11/24 | 8 / 0 | 12 / 0 | 12 / 1 |
| infant18mo | opm_dense | focal 320nAm | 23/24 | 5 / 4 | 3 / 3 | 2 / 2 |
| infant18mo | opm_dense | patch 80nAm | 8/24 | 9 / 0 | 14 / 0 | 9 / 0 |
| infant18mo | opm_dense | patch 320nAm | 23/24 | 3 / 2 | 3 / 2 | 3 / 3 |
| infant18mo | opm_matched | focal 80nAm | 11/24 | 13 / 2 | 12 / 1 | 14 / 2 |
| infant18mo | opm_matched | focal 320nAm | 23/24 | 4 / 3 | 3 / 2 | 3 / 2 |
| infant18mo | opm_matched | patch 80nAm | 7/24 | 13 / 0 | 14 / 2 | 15 / 1 |
| infant18mo | opm_matched | patch 320nAm | 22/24 | 2 / 1 | 4 / 2 | 5 / 4 |
| infant18mo | squid | focal 80nAm | 9/24 | 13 / 1 | 13 / 0 | 15 / 1 |
| infant18mo | squid | focal 320nAm | 23/24 | 3 / 2 | 7 / 7 | 5 / 4 |
| infant18mo | squid | patch 80nAm | 8/24 | 13 / 1 | 19 / 4 | 16 / 2 |
| infant18mo | squid | patch 320nAm | 21/24 | 2 / 1 | 4 / 2 | 7 / 4 |
| infant12mo | opm_dense | focal 80nAm | 13/24 | 5 / 0 | 7 / 0 | 9 / 0 |
| infant12mo | opm_dense | focal 320nAm | 24/24 | 0 / 0 | 1 / 1 | 1 / 1 |
| infant12mo | opm_dense | patch 80nAm | 11/24 | 10 / 1 | 11 / 1 | 7 / 0 |
| infant12mo | opm_dense | patch 320nAm | 23/24 | 1 / 1 | 1 / 1 | 2 / 1 |
| infant12mo | opm_matched | focal 80nAm | 13/24 | 8 / 0 | 6 / 1 | 10 / 1 |
| infant12mo | opm_matched | focal 320nAm | 23/24 | 2 / 2 | 1 / 1 | 1 / 0 |
| infant12mo | opm_matched | patch 80nAm | 9/24 | 9 / 0 | 7 / 1 | 13 / 1 |
| infant12mo | opm_matched | patch 320nAm | 23/24 | 1 / 1 | 0 / 0 | 3 / 2 |
| infant12mo | squid | focal 80nAm | 12/24 | 9 / 0 | 17 / 5 | 14 / 4 |
| infant12mo | squid | focal 320nAm | 23/24 | 0 / 0 | 0 / 0 | 3 / 2 |
| infant12mo | squid | patch 80nAm | 8/24 | 9 / 1 | 18 / 6 | 14 / 2 |
| infant12mo | squid | patch 320nAm | 23/24 | 1 / 1 | 5 / 4 | 3 / 2 |

Totals over the six anatomies (576 events per array): gross ECD errors 205 (Neuromag), 199 (matched OPM), 177 (dense OPM), of which among detected events 24, 24, 24; gross dSPM errors 262, 203, 200 (detected: 63, 28, 29). Of all gross errors (both estimators, three arrays), 85 % belong to undetected events and 76 % to 80-nAm sources, whose estimates are largely those of noise.
