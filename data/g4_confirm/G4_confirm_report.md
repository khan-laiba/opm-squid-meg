# G4 confirmatory spike-detection run

Status: confirmatory (all nine anatomies, one clean commit, declared settings); commits ae458a9.

Endpoint (fixed before the run, configs/g4_confirmatory.toml): opm_dense/opm vs squid/combined, practical detector at 1 false event per minute (thresholds frozen on the calibration null), focal spikes 10-20 mm, noise replicate 0; two-sided exact sign-flip test on the per-location differences in detection counts (detection.sign_flip_p); Holm over the nine anatomies (alpha 0.05); effect: paired S50 ratio Neuromag / OPM from the location-pooled detection curves, location bootstrap (censored interval). New seeds, newly drawn locations (vertex-disjoint from the exploratory ones), independent null data for fitting thresholds and for evaluating false-event rates.

## Endpoint

| anatomy | S50 Neuromag / dense [nAm] | S50 ratio [95 % CI] | locations OPM / Neuromag | p | Holm p | exploratory: ratio [CI], locations, p |
|---|---|---|---|---|---|---|
| adult | 50 / 35 | 1.45 [1.21-1.68] | 24 / 0 | 0 | 0 | 1.36 [1.10-1.58], 11/1, p 0.00391 |
| school | 45 / 31 | 1.44 [1.26-1.58] | 25 / 1 | 0 | 0 | 1.60 [1.31-1.87], 13/0, p 0.000244 |
| size2yr | 54 / 37 | 1.47 [1.24-1.62] | 27 / 1 | 0 | 0 | 1.65 [1.41-1.87], 17/0, p 1.53e-05 |
| infant2yr | 35 / 27 | 1.28 [1.18-1.45] | 28 / 2 | 0 | 0 | 1.44 [1.20-1.65], 14/2, p 0.00153 |
| infant18mo | 41 / 29 | 1.43 [1.26-1.61] | 30 / 1 | 0 | 0 | 1.33 [1.16-1.62], 16/0, p 3.05e-05 |
| infant12mo | 41 / 29 | 1.39 [1.25-1.56] | 29 / 1 | 0 | 0 | 1.46 [1.21-1.84], 14/1, p 0.000488 |
| childA | 51 / 40 | 1.27 [1.07-1.38] | 14 / 3 | 0.0103 | 0.0103 | 1.27 [0.93-1.63], 11/4, p 0.0225 |
| childB | 57 / 42 | 1.35 [1.15-1.56] | 21 / 2 | 0 | 0 | 1.18 [1.02-1.51], 12/3, p 0.00854 |
| childC | 50 / 39 | 1.30 [1.15-1.43] | 25 / 2 | 0 | 0 | 1.43 [1.06-1.62], 14/2, p 0.00372 |

Holm over the anatomies: 9 of 9 below 0.05.

## Monte Carlo variability (independent noise realizations of the same events)

Anatomies passing Holm per replicate: [9, 9, 9, 9, 9] (all pass in 5 of 5 replicates).

| anatomy | S50 ratio per replicate | p per replicate | SD log2 ratio | pooled ratio [CI] |
|---|---|---|---|---|
| adult | 1.45, 1.46, 1.36, 1.39, 1.41 | 0, 0, 0, 1.9e-06, 0 | 0.041 | 1.41 [1.22-1.56] |
| school | 1.44, 1.36, 1.49, 1.49, 1.39 | 0, 0, 0, 0, 0 | 0.060 | 1.44 [1.32-1.53] |
| size2yr | 1.47, 1.38, 1.44, 1.43, 1.34 | 0, 0, 0, 0, 0 | 0.052 | 1.42 [1.25-1.53] |
| infant2yr | 1.28, 1.21, 1.23, 1.21, 1.37 | 0, 0, 0, 0, 0 | 0.076 | 1.25 [1.18-1.37] |
| infant18mo | 1.43, 1.30, 1.32, 1.49, 1.46 | 0, 0, 0, 0, 0 | 0.087 | 1.39 [1.30-1.53] |
| infant12mo | 1.39, 1.27, 1.31, 1.32, 1.30 | 0, 0, 0, 5e-05, 0 | 0.051 | 1.31 [1.22-1.46] |
| childA | 1.27, 1.16, 1.24, 1.26, 1.19 | 0.01, 0.031, 0.026, 0.0059, 0.011 | 0.057 | 1.22 [1.09-1.33] |
| childB | 1.35, 1.33, 1.25, 1.15, 1.29 | 0, 0, 0.00035, 0.0015, 0.001 | 0.089 | 1.27 [1.12-1.44] |
| childC | 1.30, 1.24, 1.27, 1.35, 1.25 | 0, 0.0001, 0.00015, 0, 0 | 0.049 | 1.28 [1.17-1.38] |

## Secondary families (Holm over the anatomies within each family; not confirmatory)

| family | anatomies passing | largest Holm p |
|---|---|---|
| opm_dense/opm_vs_squid/combined/practical@1/pooled | 9 / 9 | 0.0015 |
| opm_dense/opm_vs_squid/combined/practical@1_matched/replicate0 | 9 / 9 | 0.0103 |
| opm_dense/opm_vs_squid/combined/practical@1_matched/pooled | 9 / 9 | 0.0024 |
| opm_dense/opm_vs_squid/combined/oracle/replicate0 | 9 / 9 | 0.00015 |
| opm_dense/opm_vs_squid/combined/oracle/pooled | 9 / 9 | 0 |
| opm_dense/opm_vs_squid/combined/mismatch@1/replicate0 | 9 / 9 | 0.026 |
| opm_dense/opm_vs_squid/combined/mismatch@1/pooled | 9 / 9 | 0.00255 |
| opm_dense/opm_vs_squid/combined/mismatch@1_matched/replicate0 | 8 / 9 | 0.0704 |
| opm_dense/opm_vs_squid/combined/mismatch@1_matched/pooled | 9 / 9 | 0.00115 |
| opm_matched/opm_vs_squid/combined/practical@1/replicate0 | 5 / 9 | 1 |
| opm_matched/opm_vs_squid/combined/practical@1/pooled | 7 / 9 | 0.391 |
| opm_matched/opm_vs_squid/combined/practical@1_matched/replicate0 | 5 / 9 | 1 |
| opm_matched/opm_vs_squid/combined/practical@1_matched/pooled | 7 / 9 | 0.878 |
| opm_matched/opm_vs_squid/combined/oracle/replicate0 | 4 / 9 | 1 |
| opm_matched/opm_vs_squid/combined/oracle/pooled | 7 / 9 | 1 |
| opm_matched/opm_vs_squid/combined/mismatch@1/replicate0 | 6 / 9 | 1 |
| opm_matched/opm_vs_squid/combined/mismatch@1/pooled | 6 / 9 | 0.929 |
| opm_matched/opm_vs_squid/combined/mismatch@1_matched/replicate0 | 6 / 9 | 1 |
| opm_matched/opm_vs_squid/combined/mismatch@1_matched/pooled | 6 / 9 | 0.732 |

## Detector mismatch (templates at the geometric midpoints of the injected stretches; dictionary with a 1-layer BEM and a 2-mm / 2-deg coregistration error)

| anatomy | ratio, mismatch [CI] | locations, p | ratio change (Neuromag cost / OPM cost) [CI] | cost Neuromag | cost dense |
|---|---|---|---|---|---|
| adult | 1.50 [1.30-1.67] | 27/0, p 0 | 1.03 [0.91-1.16] | 1.09 [1.01-1.22] | 1.05 [0.99-1.17] |
| school | 1.39 [1.23-1.52] | 28/2, p 0 | 0.97 [0.89-1.05] | 1.03 [0.97-1.11] | 1.07 [0.99-1.14] |
| size2yr | 1.49 [1.26-1.67] | 29/0, p 0 | 1.01 [0.93-1.09] | 1.06 [1.02-1.14] | 1.05 [0.99-1.13] |
| infant2yr | 1.42 [1.26-1.68] | 29/0, p 0 | 1.11 [0.99-1.27] | 1.14 [1.05-1.29] | 1.03 [0.98-1.09] |
| infant18mo | 1.44 [1.25-1.65] | 28/1, p 0 | 1.01 [0.93-1.09] | 1.03 [0.98-1.09] | 1.02 [0.96-1.10] |
| infant12mo | 1.43 [1.25-1.60] | 28/2, p 0 | 1.02 [0.95-1.10] | 1.04 [0.97-1.10] | 1.01 [0.97-1.06] |
| childA | 1.23 [1.06-1.41] | 15/4, p 0.026 | 0.97 [0.89-1.08] | 1.06 [1.01-1.12] | 1.10 [1.01-1.17] |
| childB | 1.23 [1.02-1.48] | 20/7, p 0.00365 | 0.91 [0.79-1.06] | 1.02 [0.93-1.11] | 1.12 [1.00-1.25] |
| childC | 1.34 [1.18-1.50] | 26/4, p 0 | 1.03 [0.94-1.16] | 1.08 [1.02-1.15] | 1.05 [0.94-1.15] |

## False events per minute (thresholds for 1 per minute)

Frozen thresholds (calibration null) on the held-out and on the evaluation null; matched thresholds (fitted on the held-out null) on the evaluation null, which neither threshold saw.

| anatomy | array | detector | held-out, frozen | evaluation, frozen | evaluation, matched |
|---|---|---|---|---|---|
| adult | squid/combined | primary | 1.20 | 1.05 | 0.95 |
| adult | squid/combined | mismatch | 0.95 | 0.85 | 0.95 |
| adult | opm_matched/opm | primary | 1.20 | 1.30 | 1.05 |
| adult | opm_matched/opm | mismatch | 1.65 | 0.80 | 0.50 |
| adult | opm_dense/opm | primary | 1.10 | 1.20 | 1.20 |
| adult | opm_dense/opm | mismatch | 1.50 | 1.40 | 0.95 |
| school | squid/combined | primary | 1.60 | 1.35 | 0.70 |
| school | squid/combined | mismatch | 1.30 | 1.20 | 1.15 |
| school | opm_matched/opm | primary | 1.00 | 0.65 | 0.70 |
| school | opm_matched/opm | mismatch | 1.20 | 1.00 | 1.00 |
| school | opm_dense/opm | primary | 0.90 | 1.25 | 1.60 |
| school | opm_dense/opm | mismatch | 0.95 | 0.80 | 1.10 |
| size2yr | squid/combined | primary | 1.15 | 1.45 | 1.25 |
| size2yr | squid/combined | mismatch | 1.20 | 1.00 | 0.85 |
| size2yr | opm_matched/opm | primary | 1.30 | 1.55 | 1.15 |
| size2yr | opm_matched/opm | mismatch | 1.05 | 0.80 | 0.80 |
| size2yr | opm_dense/opm | primary | 1.40 | 1.40 | 1.00 |
| size2yr | opm_dense/opm | mismatch | 1.00 | 1.00 | 1.00 |
| infant2yr | squid/combined | primary | 1.75 | 2.25 | 1.45 |
| infant2yr | squid/combined | mismatch | 1.00 | 1.05 | 1.20 |
| infant2yr | opm_matched/opm | primary | 1.10 | 1.15 | 0.80 |
| infant2yr | opm_matched/opm | mismatch | 1.30 | 1.25 | 1.15 |
| infant2yr | opm_dense/opm | primary | 1.20 | 1.25 | 1.15 |
| infant2yr | opm_dense/opm | mismatch | 1.30 | 1.20 | 0.85 |
| infant18mo | squid/combined | primary | 0.75 | 0.80 | 1.15 |
| infant18mo | squid/combined | mismatch | 0.65 | 0.75 | 1.15 |
| infant18mo | opm_matched/opm | primary | 1.25 | 1.00 | 0.70 |
| infant18mo | opm_matched/opm | mismatch | 1.75 | 1.45 | 0.90 |
| infant18mo | opm_dense/opm | primary | 0.65 | 1.15 | 1.25 |
| infant18mo | opm_dense/opm | mismatch | 1.35 | 1.65 | 1.10 |
| infant12mo | squid/combined | primary | 1.35 | 1.05 | 0.90 |
| infant12mo | squid/combined | mismatch | 1.40 | 1.30 | 1.05 |
| infant12mo | opm_matched/opm | primary | 1.15 | 0.75 | 0.45 |
| infant12mo | opm_matched/opm | mismatch | 1.50 | 1.30 | 0.95 |
| infant12mo | opm_dense/opm | primary | 0.75 | 1.00 | 1.40 |
| infant12mo | opm_dense/opm | mismatch | 0.75 | 1.10 | 1.45 |
| childA | squid/combined | primary | 0.95 | 0.75 | 0.80 |
| childA | squid/combined | mismatch | 1.50 | 1.65 | 0.95 |
| childA | opm_matched/opm | primary | 1.40 | 1.35 | 1.05 |
| childA | opm_matched/opm | mismatch | 1.40 | 1.60 | 1.25 |
| childA | opm_dense/opm | primary | 1.15 | 0.90 | 0.90 |
| childA | opm_dense/opm | mismatch | 1.00 | 1.00 | 1.10 |
| childB | squid/combined | primary | 0.70 | 0.70 | 1.35 |
| childB | squid/combined | mismatch | 0.90 | 0.80 | 1.00 |
| childB | opm_matched/opm | primary | 1.65 | 1.55 | 1.20 |
| childB | opm_matched/opm | mismatch | 0.95 | 0.85 | 0.90 |
| childB | opm_dense/opm | primary | 1.55 | 1.40 | 1.05 |
| childB | opm_dense/opm | mismatch | 2.05 | 1.65 | 0.90 |
| childC | squid/combined | primary | 1.10 | 1.40 | 1.35 |
| childC | squid/combined | mismatch | 1.40 | 1.30 | 0.90 |
| childC | opm_matched/opm | primary | 0.75 | 1.65 | 2.05 |
| childC | opm_matched/opm | mismatch | 1.60 | 1.50 | 1.15 |
| childC | opm_dense/opm | primary | 1.60 | 2.15 | 1.50 |
| childC | opm_dense/opm | mismatch | 1.40 | 1.20 | 0.80 |

Simulated spikes in a model; the detector knows the morphology family and (primary variant) the forward model; no claim about clinical detection is made from this run alone.
