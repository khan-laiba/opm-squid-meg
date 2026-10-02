# G4: practical detector at matched held-out false-event rates (check)

The frozen thresholds (set on 20 min of calibration null data) give unequal false-event rates on the 20 min of held-out null data. Here every detector's threshold is set on the held-out null to 1 false event per minute and the stored simulations are re-evaluated (no new simulation; in-sample for the held-out data). Paired dense OPM vs Neuromag combined and matched OPM vs Neuromag combined: locations favouring OPM / Neuromag, exact sign-flip p (uncorrected), strength ratio Neuromag / OPM [95 % location bootstrap].

| anatomy | pair | thresholds | 10-20 mm | 20-30 mm | 30-45 mm | 45-70 mm |
|---|---|---|---|---|---|---|
| adult | opm_dense vs squid/combined | frozen | 8/2, p 0.037, 1.29 [1.03-1.51] | 5/1, p 0.19, 1.07 [0.95-1.17] | 4/0, p 0.12, 1.01 [0.93-1.08] | 0/2, p 0.5, 0.98 [0.93-1.00] |
| adult | opm_dense vs squid/combined | matched | 9/1, p 0.012, 1.32 [1.04-1.55] | 7/1, p 0.062, 1.10 [0.98-1.22] | 5/0, p 0.062, 1.03 [0.95-1.10] | 0/2, p 0.5, 0.98 [0.93-1.00] |
| adult | opm_matched vs squid/combined | frozen | 8/4, p 0.52, 1.08 [0.84-1.34] | 4/6, p 0.81, 0.98 [0.85-1.12] | 1/3, p 0.62, 0.94 [0.83-1.01] | 1/9, p 0.016, 0.82 [open-0.93] |
| adult | opm_matched vs squid/combined | matched | 8/4, p 0.44, 1.11 [0.84-1.39] | 5/5, p 1, 1.02 [0.86-1.14] | 1/2, p 1, 0.96 [0.86-1.02] | 1/8, p 0.027, 0.84 [open-0.95] |
| school | opm_dense vs squid/combined | frozen | 11/1, p 0.0029, 1.40 [1.16-1.60] | 6/2, p 0.34, 1.07 [0.97-1.17] | 6/3, p 0.4, 1.10 [0.96-1.26] | 1/4, p 0.38, 0.96 [open-1.04] |
| school | opm_dense vs squid/combined | matched | 12/1, p 0.0012, 1.42 [1.16-1.64] | 7/1, p 0.055, 1.09 [1.02-1.18] | 7/3, p 0.27, 1.13 [0.97-1.29] | 4/2, p 0.69, 1.04 [0.98-open] |
| school | opm_matched vs squid/combined | frozen | 6/5, p 0.63, 1.15 [0.85-1.46] | 4/3, p 0.83, 1.03 [0.95-1.15] | 2/4, p 0.41, 0.90 [0.78-1.04] | 0/7, p 0.016, 0.87 [open-0.96] |
| school | opm_matched vs squid/combined | matched | 6/4, p 0.39, 1.16 [0.90-1.46] | 4/3, p 0.83, 1.03 [0.95-1.15] | 2/4, p 0.53, 0.93 [0.83-1.04] | 0/7, p 0.016, 0.87 [open-0.96] |
| size2yr | opm_dense vs squid/combined | frozen | 16/0, p 3.1e-05, 1.61 [1.33-1.85] | 9/2, p 0.045, 1.19 [0.98-1.46] | 9/2, p 0.045, 1.13 [0.99-1.35] | 9/1, p 0.02, 1.12 [1.05-open] |
| size2yr | opm_dense vs squid/combined | matched | 15/0, p 6.1e-05, 1.59 [1.32-1.84] | 9/2, p 0.061, 1.16 [0.96-1.40] | 7/2, p 0.11, 1.09 [0.96-1.29] | 8/1, p 0.035, 1.12 [1.04-open] |
| size2yr | opm_matched vs squid/combined | frozen | 12/0, p 0.00049, 1.55 [1.11-1.78] | 7/1, p 0.055, 1.09 [1.00-1.26] | 4/2, p 0.53, 1.07 [0.98-1.20] | 3/4, p 0.47, 0.98 [open-open] |
| size2yr | opm_matched vs squid/combined | matched | 12/0, p 0.00049, 1.54 [1.13-1.76] | 7/1, p 0.055, 1.09 [0.99-1.22] | 3/2, p 1, 1.02 [0.97-1.09] | 2/3, p 0.56, 0.98 [open-open] |
| infant2yr | opm_dense vs squid/combined | frozen | 14/2, p 0.0015, 1.44 [1.20-1.65] | 6/1, p 0.22, 1.07 [0.97-1.23] | 7/3, p 0.46, 1.03 [0.93-1.15] | 7/2, p 0.15, 1.05 [0.96-1.13] |
| infant2yr | opm_dense vs squid/combined | matched | 13/2, p 0.0026, 1.40 [1.18-1.62] | 6/3, p 0.49, 1.05 [0.94-1.21] | 4/3, p 1, 0.98 [0.91-1.07] | 7/4, p 0.43, 1.03 [0.93-1.12] |
| infant2yr | opm_matched vs squid/combined | frozen | 10/5, p 0.15, 1.22 [1.00-1.42] | 8/2, p 0.27, 1.01 [0.88-1.13] | 3/5, p 0.44, 0.95 [0.82-1.06] | 5/4, p 1, 1.01 [0.93-1.08] |
| infant2yr | opm_matched vs squid/combined | matched | 9/7, p 0.66, 1.16 [0.92-1.39] | 7/3, p 0.59, 0.99 [0.85-1.11] | 2/6, p 0.19, 0.93 [0.80-1.04] | 3/5, p 0.44, 0.97 [0.86-1.04] |
| infant18mo | opm_dense vs squid/combined | frozen | 16/0, p 3.1e-05, 1.33 [1.16-1.62] | 6/2, p 0.19, 1.05 [0.89-1.25] | 8/3, p 0.18, 1.05 [0.98-1.13] | 7/1, p 0.047, Neuromag does not reach 50 % [1.05-open] |
| infant18mo | opm_dense vs squid/combined | matched | 16/0, p 3.1e-05, 1.33 [1.16-1.62] | 6/2, p 0.19, 1.05 [0.89-1.25] | 7/3, p 0.27, 1.04 [0.97-1.12] | 7/1, p 0.15, Neuromag does not reach 50 % [1.05-open] |
| infant18mo | opm_matched vs squid/combined | frozen | 12/0, p 0.00049, 1.10 [1.02-1.28] | 9/6, p 0.45, 1.03 [0.78-1.33] | 4/3, p 0.77, 1.01 [0.95-1.11] | 2/4, p 0.53, neither reaches 50 % [open-open] |
| infant18mo | opm_matched vs squid/combined | matched | 11/1, p 0.0044, 1.04 [0.94-1.20] | 9/6, p 0.59, 1.03 [0.73-1.33] | 3/4, p 1, 0.98 [0.93-1.05] | 2/4, p 0.53, neither reaches 50 % [open-open] |
| infant12mo | opm_dense vs squid/combined | frozen | 14/1, p 0.00049, 1.46 [1.21-1.84] | 9/3, p 0.22, 1.12 [1.01-1.25] | 8/1, p 0.031, 1.10 [1.03-1.23] | 3/2, p 0.75, 1.02 [0.94-1.15] |
| infant12mo | opm_dense vs squid/combined | matched | 14/1, p 0.00043, 1.51 [1.24-1.87] | 11/3, p 0.11, 1.15 [1.03-1.29] | 9/1, p 0.018, 1.11 [1.04-1.24] | 3/2, p 0.75, 1.02 [0.94-1.15] |
| infant12mo | opm_matched vs squid/combined | frozen | 11/4, p 0.051, 1.28 [1.09-1.46] | 9/3, p 0.097, 1.12 [1.01-1.25] | 4/4, p 0.8, 1.01 [0.93-1.11] | 2/2, p 1, 0.99 [0.93-1.06] |
| infant12mo | opm_matched vs squid/combined | matched | 11/4, p 0.033, 1.34 [1.12-1.53] | 9/3, p 0.078, 1.13 [1.02-1.25] | 6/3, p 0.31, 1.05 [0.97-1.17] | 2/2, p 1, 0.99 [0.93-1.06] |

Held-out false events per minute at the frozen thresholds (matched: 1.00 by construction):

| anatomy | squid/combined | opm_matched/opm | opm_dense/opm |
|---|---|---|---|
| adult | 1.00 | 0.65 | 0.80 |
| school | 1.00 | 0.70 | 0.40 |
| size2yr | 0.45 | 0.70 | 1.25 |
| infant2yr | 0.70 | 1.55 | 1.30 |
| infant18mo | 0.90 | 1.45 | 1.10 |
| infant12mo | 1.25 | 0.95 | 0.85 |
