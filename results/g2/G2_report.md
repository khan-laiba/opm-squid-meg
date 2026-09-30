# G2 report: realistic adult OPM vs Neuromag comparison (NEW)

Generated from `results/g2/g2_summary.json` (code commit 29bea4a, MNE 1.13.2). Methods: `docs/methods.md` section 8; assumptions in `docs/provenance_register.md`. This is a proposed study; no author of the reproduced papers has reviewed it.

## Common setup

- Anatomy: MNE sample subject, measured head position; 8124 target dipoles (10 nAm, cortical normal, usable oct-6 vertices); background grid of 1800 area-weighted sources.
- Band 1-40 Hz (ENBW 35.1 Hz). Intrinsic noise: SQUID magnetometers 20.7 fT, gradiometers 21.3 fT/cm; OPM at 15 fT/sqrt(Hz) 89 fT (RMS in band).
- Brain noise calibrated on good gradiometers to 37.1 fT/cm (task baseline minus empty room); predicted magnetometer level 193 fT vs 262 fT measured.
- Room field: explains 94 % (magnetometers) and 1.6 % (gradiometers) of the empty-room variance; model vs measured empty room 115 vs 115 fT, 21.4 vs 20.2 fT/cm.
- Conditions: intrinsic, intrinsic+brain, intrinsic+brain+env, projected. Detectability = known-topography matched-filter SNR sqrt(s^T C^+ s); not an event detection rate or a localization accuracy.

## Neuromag magnetometers (102)

- Geometry: 102 channels, 1 mag + 2 planar grad per site; scalp-to-sensor distance median 29.8 mm (5-95 %: 25.4-35.9 mm).
- Retained rank of the full array: intrinsic 306, intrinsic+brain 306, intrinsic+brain+env 306, projected 298.
- Noise composition (median RMS): intrinsic 20.7 fT, brain 192.7 fT, room 113.0 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.27 | 1.28 | 0.98 | 0.67 | 0.46 | 0.33 | 0.26 | 0.22 | 0.16 | 0.15 | 0.12 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.24 | 1.26 | 0.95 | 0.65 | 0.44 | 0.31 | 0.24 | 0.19 | 0.15 | 0.12 | 0.10 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.95 (10 s), 0.99 (60 s), projected: 0.93 (10 s), 0.98 (60 s).

## Neuromag planar gradiometers (204)

- Geometry: 204 channels, 1 mag + 2 planar grad per site; scalp-to-sensor distance median 29.8 mm (5-95 %: 25.4-35.9 mm).
- Retained rank of the full array: intrinsic 306, intrinsic+brain 306, intrinsic+brain+env 306, projected 298.
- Noise composition (median RMS): intrinsic 21.3 fT/cm, brain 37.1 fT/cm, room 1.7 fT/cm.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.17 | 1.20 | 0.91 | 0.64 | 0.43 | 0.31 | 0.25 | 0.20 | 0.15 | 0.13 | 0.11 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.17 | 1.20 | 0.91 | 0.63 | 0.43 | 0.30 | 0.24 | 0.20 | 0.15 | 0.12 | 0.11 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.89 (10 s), 0.98 (60 s), projected: 0.89 (10 s), 0.98 (60 s).

## Neuromag combined (306)

- Geometry: 306 channels, 1 mag + 2 planar grad per site; scalp-to-sensor distance median 29.8 mm (5-95 %: 25.4-35.9 mm).
- Retained rank of the full array: intrinsic 306, intrinsic+brain 306, intrinsic+brain+env 306, projected 298.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.40 | 1.38 | 1.02 | 0.69 | 0.47 | 0.34 | 0.27 | 0.23 | 0.17 | 0.15 | 0.13 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.39 | 1.37 | 1.01 | 0.68 | 0.46 | 0.32 | 0.26 | 0.21 | 0.16 | 0.14 | 0.11 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.85 (10 s), 0.97 (60 s), projected: 0.87 (10 s), 0.97 (60 s).

## OPM matched sites (99)

- Geometry: 99 channels, 99 sites, 1 (scalp normal); scalp-to-sensor distance median 7.0 mm (5-95 %: 6.3-7.0 mm); role: matched-site coverage control.
- Retained rank: intrinsic 99, intrinsic+brain 99, intrinsic+brain+env 99, projected 91.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 520.2 fT, room 115.3 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.56 | 1.54 | 1.09 | 0.71 | 0.47 | 0.33 | 0.26 | 0.21 | 0.16 | 0.14 | 0.12 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.53 | 1.51 | 1.05 | 0.67 | 0.43 | 0.29 | 0.22 | 0.17 | 0.13 | 0.10 | 0.09 |

OPM / Neuromag, median detectability ratio (bootstrap 95 % CI over targets):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.54x [0.53-0.54], OPM higher for 0 % | 2.29x [2.28-2.30], OPM higher for 100 % | 0.55x [0.55-0.56], OPM higher for 0 % |
| intrinsic+brain | 1.01x [1.01-1.02], OPM higher for 59 % | 1.14x [1.14-1.14], OPM higher for 92 % | 1.05x [1.04-1.05], OPM higher for 75 % |
| intrinsic+brain+env | 1.01x [1.01-1.02], OPM higher for 58 % | 1.12x [1.11-1.12], OPM higher for 88 % | 1.08x [1.08-1.08], OPM higher for 89 % |
| projected | 0.96x [0.96-0.97], OPM higher for 39 % | 1.05x [1.04-1.06], OPM higher for 61 % | 1.02x [1.02-1.03], OPM higher for 58 % |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.95 (10 s), 0.99 (60 s), projected: 0.95 (10 s), 0.99 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.02x, 10 mm 1.02x, 20 mm 1.02x.
- Sensitivity (vs combined, intrinsic+brain): opm asd 7fT 1.09x; opm asd 10fT 1.05x; opm asd 15fT 1.01x; opm asd 20fT 0.99x; opm asd 30fT 0.93x; background corr 5mm 1.02x; background corr 10mm 1.03x; background mag calibrated 1.01x; head x+5mm 1.01x; head x-5mm 1.01x; head y+5mm 1.02x; head y-5mm 1.00x; head z+5mm 1.00x; head z-5mm 1.03x; head pitch+5deg 1.01x; head pitch-5deg 1.02x; head well fitted 0.99x; gap 3mm 0.99x; gap 6mm 0.96x.
- By lobe (vs combined, intrinsic+brain): frontal 1.04x, parietal 1.04x, temporal 1.00x, occipital 1.02x, cingulate 0.99x, insula 0.96x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.01x, 1-10Hz 1.00x, 8-30Hz 1.02x, 30-80Hz 1.00x.

## OPM channel-budget control (204)

- Geometry: 204 channels, 204 sites, 1 (scalp normal); scalp-to-sensor distance median 7.0 mm (5-95 %: 6.4-7.0 mm); role: channel-budget control vs 204 gradiometers.
- Retained rank: intrinsic 204, intrinsic+brain 204, intrinsic+brain+env 204, projected 196.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 501.9 fT, room 114.2 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.19 | 1.88 | 1.25 | 0.78 | 0.51 | 0.37 | 0.29 | 0.24 | 0.19 | 0.17 | 0.14 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.18 | 1.87 | 1.24 | 0.77 | 0.49 | 0.35 | 0.27 | 0.22 | 0.17 | 0.14 | 0.12 |

OPM / Neuromag, median detectability ratio (bootstrap 95 % CI over targets):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.73x [0.73-0.73], OPM higher for 3 % | 3.09x [3.08-3.10], OPM higher for 100 % | 0.75x [0.75-0.76], OPM higher for 5 % |
| intrinsic+brain | 1.16x [1.16-1.16], OPM higher for 100 % | 1.35x [1.34-1.35], OPM higher for 100 % | 1.21x [1.21-1.22], OPM higher for 100 % |
| intrinsic+brain+env | 1.17x [1.17-1.18], OPM higher for 100 % | 1.32x [1.32-1.33], OPM higher for 100 % | 1.27x [1.26-1.27], OPM higher for 100 % |
| projected | 1.17x [1.16-1.17], OPM higher for 98 % | 1.31x [1.30-1.31], OPM higher for 100 % | 1.26x [1.26-1.27], OPM higher for 100 % |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.91 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.16x, 10 mm 1.16x, 20 mm 1.15x.
- Sensitivity (vs combined, intrinsic+brain): opm asd 7fT 1.33x; opm asd 10fT 1.25x; opm asd 15fT 1.16x; opm asd 20fT 1.10x; opm asd 30fT 1.01x; background corr 5mm 1.17x; background corr 10mm 1.19x; background mag calibrated 1.17x; head x+5mm 1.15x; head x-5mm 1.16x; head y+5mm 1.17x; head y-5mm 1.14x; head z+5mm 1.13x; head z-5mm 1.19x; head pitch+5deg 1.15x; head pitch-5deg 1.17x; head well fitted 1.12x; gap 3mm 1.09x; gap 6mm 1.04x.
- By lobe (vs combined, intrinsic+brain): frontal 1.22x, parietal 1.16x, temporal 1.15x, occipital 1.11x, cingulate 1.09x, insula 1.10x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.16x, 1-10Hz 1.17x, 8-30Hz 1.15x, 30-80Hz 1.11x.

## OPM densest feasible array (216)

- Geometry: 216 channels, 216 sites, 1 (scalp normal); scalp-to-sensor distance median 7.0 mm (5-95 %: 6.4-7.0 mm); role: densest feasible single-axis array (full system).
- Retained rank: intrinsic 216, intrinsic+brain 216, intrinsic+brain+env 216, projected 208.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 500.6 fT, room 114.5 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.22 | 1.91 | 1.26 | 0.78 | 0.51 | 0.37 | 0.29 | 0.24 | 0.19 | 0.17 | 0.14 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.20 | 1.90 | 1.25 | 0.77 | 0.50 | 0.35 | 0.27 | 0.22 | 0.17 | 0.14 | 0.12 |

OPM / Neuromag, median detectability ratio (bootstrap 95 % CI over targets):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.74x [0.74-0.75], OPM higher for 5 % | 3.16x [3.15-3.16], OPM higher for 100 % | 0.77x [0.77-0.77], OPM higher for 7 % |
| intrinsic+brain | 1.17x [1.17-1.17], OPM higher for 100 % | 1.36x [1.35-1.37], OPM higher for 100 % | 1.23x [1.22-1.23], OPM higher for 100 % |
| intrinsic+brain+env | 1.18x [1.18-1.19], OPM higher for 100 % | 1.33x [1.33-1.34], OPM higher for 100 % | 1.28x [1.27-1.28], OPM higher for 100 % |
| projected | 1.18x [1.17-1.18], OPM higher for 98 % | 1.32x [1.31-1.32], OPM higher for 100 % | 1.27x [1.27-1.28], OPM higher for 100 % |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.90 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.17x, 10 mm 1.17x, 20 mm 1.15x.
- Sensitivity (vs combined, intrinsic+brain): opm asd 7fT 1.34x; opm asd 10fT 1.27x; opm asd 15fT 1.17x; opm asd 20fT 1.10x; opm asd 30fT 1.02x; background corr 5mm 1.18x; background corr 10mm 1.20x; background mag calibrated 1.18x; head x+5mm 1.16x; head x-5mm 1.17x; head y+5mm 1.18x; head y-5mm 1.15x; head z+5mm 1.14x; head z-5mm 1.20x; head pitch+5deg 1.16x; head pitch-5deg 1.18x; head well fitted 1.13x; gap 3mm 1.10x; gap 6mm 1.04x.
- By lobe (vs combined, intrinsic+brain): frontal 1.24x, parietal 1.16x, temporal 1.16x, occipital 1.12x, cingulate 1.10x, insula 1.10x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.17x, 1-10Hz 1.18x, 8-30Hz 1.16x, 30-80Hz 1.12x.

## Convergence

- Background grid vs every usable vertex: median log2 ratios change by <= 0.014.
- BEM 5,120 vs 20,480 triangles (1 layer): <= 0.005; 3 vs 1 layer: <= 0.066.
- Oct-6 vs random full-resolution targets: <= 0.018; Neuromag 4-point vs accurate integration: 0.6 % in the gains.

## Limitations

- Bootstrap CIs resample cortical target locations of one anatomy; they do not include model or between-subject uncertainty, which the sensitivity analyses bound.
- Head-position variants move the head in the fixed Neuromag helmet; the room field is kept in head coordinates (8-term model; the change over a 5-mm move is second order).
- Detectability is a known-topography matched-filter SNR, not an event detection rate or localization accuracy.
- One adult anatomy and one measured head position; between-subject variability is not represented.
- OPM intrinsic noise is a declared sweep, not a device specification; OPM movement artefacts, cross-talk and calibration errors are not modelled.
