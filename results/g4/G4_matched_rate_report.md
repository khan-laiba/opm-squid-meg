# G4: practical detector at matched held-out false-event rates (check)

The frozen thresholds (set on 20 min of calibration null data) give unequal false-event rates on the 20 min of held-out null data. Here every detector's threshold is set on the held-out null to 1 false event per minute and the stored simulations are re-evaluated (no new simulation; in-sample for the held-out data). Paired dense OPM vs Neuromag combined and matched OPM vs Neuromag combined: locations favouring OPM / Neuromag, exact sign-flip p (uncorrected), strength ratio Neuromag / OPM [95 % location bootstrap].

| anatomy | pair | thresholds | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|---|---|
| adult | opm_dense vs squid/combined | frozen | 16/0, p 3.1e-05, 1.51 [1.25-1.66] | 6/2, p 0.23, 1.08 [1.00-1.20] | 6/1, p 0.28, 1.04 [0.99-1.12] | 6/1, p 0.12, 1.07 [1.00-1.16] |
| adult | opm_dense vs squid/combined | matched | 16/0, p 3.1e-05, 1.51 [1.25-1.66] | 6/2, p 0.23, 1.08 [1.00-1.20] | 7/1, p 0.15, 1.07 [1.00-1.17] | 6/1, p 0.12, 1.07 [1.00-1.16] |
| adult | opm_matched vs squid/combined | frozen | 9/3, p 0.04, 1.27 [1.05-1.52] | 5/4, p 0.85, 1.01 [0.84-1.10] | 4/5, p 0.62, 0.96 [0.88-1.03] | 0/7, p 0.016, n/a [0.82-0.95] |
| adult | opm_matched vs squid/combined | matched | 9/2, p 0.019, 1.31 [1.05-1.55] | 6/4, p 1, 1.06 [0.93-1.19] | 4/3, p 1, 1.02 [0.94-1.13] | 0/6, p 0.031, n/a [0.82-0.96] |
| school | opm_dense vs squid/combined | frozen | 16/0, p 3.1e-05, 1.42 [1.20-1.57] | 6/3, p 0.29, 1.06 [0.96-1.21] | 4/3, p 0.77, 1.00 [0.93-1.10] | 6/1, p 0.11, 1.12 [1.04-1.22] |
| school | opm_dense vs squid/combined | matched | 16/0, p 3.1e-05, 1.39 [1.19-1.52] | 6/3, p 0.31, 1.03 [0.93-1.18] | 4/5, p 1, 0.97 [0.89-1.09] | 6/3, p 0.4, 1.08 [0.98-1.18] |
| school | opm_matched vs squid/combined | frozen | 6/3, p 0.34, 1.13 [0.95-1.39] | 2/5, p 0.45, 0.94 [0.85-1.02] | 3/6, p 0.31, 0.97 [0.83-1.07] | 2/4, p 0.31, 1.02 [0.88-1.12] |
| school | opm_matched vs squid/combined | matched | 8/3, p 0.24, 1.18 [0.96-1.42] | 4/6, p 0.75, 0.95 [0.84-1.08] | 3/8, p 0.15, 0.94 [0.80-1.04] | 2/6, p 0.15, 0.97 [0.85-1.09] |
| size2yr | opm_dense vs squid/combined | frozen | 16/0, p 3.1e-05, 1.64 [1.40-1.85] | 7/3, p 0.19, 1.12 [0.95-1.37] | 10/1, p 0.0088, 1.18 [1.06-1.31] | 8/1, p 0.031, 1.12 [1.03-1.24] |
| size2yr | opm_dense vs squid/combined | matched | 15/0, p 6.1e-05, 1.62 [1.38-1.84] | 6/3, p 0.25, 1.10 [0.93-1.33] | 9/1, p 0.016, 1.18 [1.04-1.30] | 7/1, p 0.062, 1.09 [1.01-1.20] |
| size2yr | opm_matched vs squid/combined | frozen | 12/1, p 0.0017, 1.55 [1.11-1.78] | 7/1, p 0.055, 1.09 [1.00-1.26] | 4/1, p 0.31, 1.08 [0.99-1.21] | 4/4, p 0.64, 0.98 [0.87-1.08] |
| size2yr | opm_matched vs squid/combined | matched | 12/0, p 0.00049, 1.54 [1.13-1.76] | 7/1, p 0.062, 1.06 [0.98-1.16] | 2/0, p 0.5, 1.05 [1.00-1.16] | 3/5, p 0.44, 0.96 [0.86-1.07] |
| infant2yr | opm_dense vs squid/combined | frozen | 16/2, p 0.00047, 1.50 [1.26-1.68] | 7/1, p 0.15, 1.07 [0.97-1.23] | 6/3, p 0.62, 1.01 [0.92-1.12] | 7/3, p 0.22, 1.06 [0.97-1.15] |
| infant2yr | opm_dense vs squid/combined | matched | 15/2, p 0.00079, 1.49 [1.25-1.66] | 6/4, p 0.66, 1.03 [0.92-1.19] | 5/3, p 0.8, 0.99 [0.91-1.08] | 7/4, p 0.34, 1.05 [0.96-1.14] |
| infant2yr | opm_matched vs squid/combined | frozen | 11/5, p 0.052, 1.28 [1.05-1.45] | 8/1, p 0.15, 1.01 [0.88-1.13] | 3/6, p 0.31, 0.93 [0.81-1.05] | 5/4, p 1, 1.01 [0.94-1.07] |
| infant2yr | opm_matched vs squid/combined | matched | 11/5, p 0.3, 1.24 [1.00-1.44] | 7/3, p 0.59, 0.99 [0.85-1.11] | 3/6, p 0.31, 0.93 [0.81-1.05] | 4/4, p 0.67, 0.99 [0.90-1.05] |
| infant18mo | opm_dense vs squid/combined | frozen | 16/0, p 3.1e-05, 1.33 [1.16-1.62] | 6/1, p 0.094, 1.06 [0.90-1.26] | 8/3, p 0.18, 1.05 [0.98-1.13] | 7/1, p 0.047, n/a [1.03-1.29] |
| infant18mo | opm_dense vs squid/combined | matched | 16/0, p 3.1e-05, 1.30 [1.14-1.58] | 5/2, p 0.28, 1.04 [0.89-1.22] | 7/3, p 0.46, 1.02 [0.94-1.11] | 7/1, p 0.12, n/a [1.03-1.29] |
| infant18mo | opm_matched vs squid/combined | frozen | 12/0, p 0.00049, 1.10 [1.02-1.27] | 9/5, p 0.33, 1.04 [0.78-1.34] | 4/2, p 0.53, 1.03 [0.98-1.12] | 2/4, p 0.69, n/a [0.90-1.07] |
| infant18mo | opm_matched vs squid/combined | matched | 10/1, p 0.0078, 1.04 [0.94-1.18] | 8/6, p 0.72, 1.02 [0.73-1.29] | 3/5, p 0.73, 0.97 [0.90-1.03] | 2/4, p 0.53, n/a [0.90-1.07] |
| infant12mo | opm_dense vs squid/combined | frozen | 14/1, p 0.00043, 1.39 [1.19-1.68] | 5/3, p 0.56, 1.04 [0.96-1.17] | 7/4, p 0.43, 1.07 [0.96-1.20] | 6/3, p 0.31, 1.08 [0.96-1.20] |
| infant12mo | opm_dense vs squid/combined | matched | 15/1, p 0.00024, 1.49 [1.24-1.82] | 8/2, p 0.092, 1.05 [0.99-1.19] | 7/4, p 0.43, 1.07 [0.96-1.20] | 6/3, p 0.31, 1.08 [0.96-1.20] |
| infant12mo | opm_matched vs squid/combined | frozen | 9/1, p 0.014, 1.23 [1.07-1.48] | 4/3, p 0.59, 1.03 [0.95-1.16] | 4/2, p 0.69, 1.03 [0.95-1.10] | 6/4, p 0.66, 1.05 [0.93-1.16] |
| infant12mo | opm_matched vs squid/combined | matched | 10/1, p 0.0068, 1.27 [1.10-1.50] | 5/3, p 0.44, 1.03 [0.96-1.16] | 4/2, p 0.69, 1.03 [0.95-1.10] | 6/4, p 0.66, 1.05 [0.93-1.16] |

Held-out false events per minute at the frozen thresholds (matched: 1.00 by construction):

| anatomy | squid/combined | opm_matched/opm | opm_dense/opm |
|---|---|---|---|
| adult | 0.95 | 0.50 | 0.90 |
| school | 0.55 | 0.80 | 1.00 |
| size2yr | 0.40 | 0.60 | 0.90 |
| infant2yr | 0.65 | 1.55 | 1.25 |
| infant18mo | 0.95 | 1.45 | 1.05 |
| infant12mo | 1.25 | 1.05 | 0.70 |
