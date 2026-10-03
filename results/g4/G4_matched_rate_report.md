# G4: practical detector at matched held-out false-event rates (check)

The frozen thresholds (set on 20 min of calibration null data) give unequal false-event rates on the 20 min of held-out null data. Here every detector's threshold is set on the held-out null to 1 false event per minute and the stored simulations are re-evaluated (no new simulation; in-sample for the held-out data). Paired dense OPM vs Neuromag combined and matched OPM vs Neuromag combined: locations favouring OPM / Neuromag, exact sign-flip p (uncorrected), strength ratio Neuromag / OPM [95 % location bootstrap].

| anatomy | pair | thresholds | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|---|---|
| adult | opm_dense vs squid/combined | frozen | 11/1, p 0.0039, 1.36 [1.10-1.58] | 9/1, p 0.014, 1.20 [1.06-1.34] | 5/1, p 0.19, 1.09 [1.00-1.22] | 9/0, p 0.0039, 1.19 [1.00-open] |
| adult | opm_dense vs squid/combined | matched | 12/1, p 0.0022, 1.40 [1.12-1.59] | 12/0, p 0.00049, 1.25 [1.08-1.41] | 5/1, p 0.19, 1.09 [1.00-1.22] | 9/0, p 0.0039, 1.19 [1.00-open] |
| adult | opm_matched vs squid/combined | frozen | 7/4, p 0.4, 1.16 [0.88-1.44] | 4/5, p 1, 1.03 [0.87-1.15] | 2/5, p 0.77, 0.98 [0.86-1.13] | 3/6, p 0.4, 0.94 [open-open] |
| adult | opm_matched vs squid/combined | matched | 8/4, p 0.32, 1.19 [0.89-1.47] | 8/3, p 0.23, 1.07 [0.93-1.19] | 2/5, p 0.77, 0.98 [0.86-1.13] | 5/6, p 0.79, 0.94 [open-open] |
| school | opm_dense vs squid/combined | frozen | 13/0, p 0.00024, 1.60 [1.31-1.87] | 4/1, p 0.31, 1.01 [0.94-1.12] | 8/3, p 0.24, 1.04 [0.85-1.24] | 6/3, p 0.78, 1.08 [open-open] |
| school | opm_dense vs squid/combined | matched | 14/0, p 0.00012, 1.63 [1.34-1.88] | 5/1, p 0.19, 1.04 [0.94-1.15] | 8/2, p 0.16, 1.05 [0.86-1.26] | 5/3, p 1, 1.06 [open-open] |
| school | opm_matched vs squid/combined | frozen | 8/3, p 0.25, 1.27 [0.99-1.52] | 4/4, p 1, 1.01 [0.91-1.13] | 4/6, p 1, 0.95 [0.78-1.17] | 4/5, p 0.61, 1.05 [open-open] |
| school | opm_matched vs squid/combined | matched | 10/3, p 0.064, 1.35 [1.03-1.57] | 4/4, p 0.8, 1.03 [0.91-1.21] | 5/5, p 1, 0.96 [0.79-1.18] | 4/5, p 0.75, 1.05 [open-open] |
| size2yr | opm_dense vs squid/combined | frozen | 17/0, p 1.5e-05, 1.65 [1.41-1.87] | 7/3, p 0.24, 1.11 [0.94-1.35] | 8/2, p 0.049, 1.16 [1.01-1.35] | 8/1, p 0.031, 1.11 [open-open] |
| size2yr | opm_dense vs squid/combined | matched | 15/0, p 6.1e-05, 1.56 [1.31-1.80] | 6/3, p 0.39, 1.06 [0.92-1.26] | 8/3, p 0.094, 1.12 [0.97-1.29] | 7/2, p 0.12, 1.09 [open-open] |
| size2yr | opm_matched vs squid/combined | frozen | 14/0, p 0.00012, 1.60 [1.21-1.82] | 7/1, p 0.055, 1.08 [0.99-1.22] | 4/1, p 0.31, 1.07 [1.00-1.20] | 3/4, p 0.59, 0.98 [open-open] |
| size2yr | opm_matched vs squid/combined | matched | 12/0, p 0.00049, 1.52 [1.13-1.74] | 7/1, p 0.055, 1.06 [0.98-1.16] | 5/1, p 0.22, 1.06 [1.00-1.14] | 2/3, p 0.56, 0.98 [open-open] |
| infant2yr | opm_dense vs squid/combined | frozen | 14/2, p 0.0015, 1.44 [1.20-1.65] | 6/1, p 0.22, 1.07 [0.97-1.23] | 7/3, p 0.46, 1.03 [0.93-1.15] | 7/2, p 0.15, 1.05 [0.96-1.13] |
| infant2yr | opm_dense vs squid/combined | matched | 13/2, p 0.0026, 1.40 [1.18-1.62] | 6/4, p 0.66, 1.03 [0.91-1.19] | 4/3, p 1, 0.98 [0.91-1.07] | 7/4, p 0.43, 1.03 [0.93-1.12] |
| infant2yr | opm_matched vs squid/combined | frozen | 10/5, p 0.15, 1.22 [1.00-1.42] | 8/2, p 0.27, 1.01 [0.88-1.13] | 3/5, p 0.44, 0.95 [0.82-1.06] | 5/4, p 1, 1.01 [0.93-1.09] |
| infant2yr | opm_matched vs squid/combined | matched | 9/7, p 0.66, 1.16 [0.92-1.39] | 7/3, p 0.59, 0.99 [0.85-1.11] | 2/6, p 0.19, 0.93 [0.80-1.04] | 3/5, p 0.44, 0.97 [0.86-1.04] |
| infant18mo | opm_dense vs squid/combined | frozen | 16/0, p 3.1e-05, 1.33 [1.16-1.62] | 6/2, p 0.19, 1.05 [0.89-1.25] | 8/3, p 0.18, 1.05 [0.98-1.13] | 7/1, p 0.047, > 1.13 (Neuromag does not reach 50 %) [open-open] |
| infant18mo | opm_dense vs squid/combined | matched | 16/0, p 3.1e-05, 1.33 [1.16-1.62] | 6/2, p 0.19, 1.05 [0.89-1.25] | 7/3, p 0.27, 1.04 [0.97-1.12] | 7/1, p 0.15, > 1.13 (Neuromag does not reach 50 %) [open-open] |
| infant18mo | opm_matched vs squid/combined | frozen | 12/0, p 0.00049, 1.10 [1.02-1.28] | 9/6, p 0.45, 1.03 [0.78-1.33] | 4/3, p 0.77, 1.01 [0.95-1.11] | 2/4, p 0.53, neither reaches 50 % [open-open] |
| infant18mo | opm_matched vs squid/combined | matched | 11/1, p 0.0044, 1.04 [0.94-1.20] | 9/6, p 0.59, 1.03 [0.73-1.33] | 3/4, p 1, 0.98 [0.93-1.05] | 2/4, p 0.53, neither reaches 50 % [open-open] |
| infant12mo | opm_dense vs squid/combined | frozen | 14/1, p 0.00049, 1.46 [1.21-1.84] | 9/3, p 0.22, 1.12 [1.01-1.25] | 8/1, p 0.031, 1.10 [1.03-1.23] | 3/2, p 0.75, 1.02 [0.93-open] |
| infant12mo | opm_dense vs squid/combined | matched | 14/1, p 0.00043, 1.51 [1.24-1.87] | 11/3, p 0.11, 1.15 [1.03-1.29] | 9/1, p 0.018, 1.11 [1.04-1.24] | 3/2, p 0.75, 1.02 [0.93-open] |
| infant12mo | opm_matched vs squid/combined | frozen | 11/4, p 0.051, 1.28 [1.09-1.46] | 9/3, p 0.097, 1.12 [1.01-1.25] | 4/4, p 0.8, 1.01 [0.93-1.11] | 2/2, p 1, 0.99 [open-open] |
| infant12mo | opm_matched vs squid/combined | matched | 11/4, p 0.033, 1.34 [1.12-1.53] | 9/3, p 0.078, 1.13 [1.02-1.25] | 6/3, p 0.31, 1.05 [0.97-1.17] | 2/2, p 1, 0.99 [open-open] |

Held-out false events per minute at the frozen thresholds (matched: 1.00 by construction):

| anatomy | squid/combined | opm_matched/opm | opm_dense/opm |
|---|---|---|---|
| adult | 1.30 | 0.70 | 0.95 |
| school | 1.40 | 0.50 | 1.05 |
| size2yr | 0.40 | 0.65 | 1.20 |
| infant2yr | 0.70 | 1.55 | 1.30 |
| infant18mo | 0.90 | 1.45 | 1.10 |
| infant12mo | 1.20 | 0.95 | 0.85 |
