# G2 report: realistic adult OPM vs Neuromag comparison (NEW)

Generated from `results/g2/g2_summary.json` (code commit 2aa6a2a, figures replotted at 977c57b; band supplement `g2_band_sensitivity.json` at 8e4c9e9; MNE 1.13.2). Methods: `docs/methods.md` section 8; assumptions in `docs/provenance_register.md`. This is a proposed study; no author of the reproduced papers has reviewed it.

## Common setup

- Anatomy: MNE sample subject, measured head position; 7661 target dipoles (10 nAm, cortical normal, usable oct-6 vertices); background grid of 1755 area-weighted sources.
- Band 1-40 Hz (ENBW 35.1 Hz). Intrinsic noise: SQUID magnetometers 20.7 fT, gradiometers 21.3 fT/cm; OPM at 15 fT/sqrt(Hz) 89 fT (RMS in band).
- Brain noise calibrated on good gradiometers to 37.1 fT/cm (task baseline minus empty room); predicted magnetometer level 192 fT vs 262 fT measured.
- Room field: explains 94 % (magnetometers) and 1.6 % (gradiometers) of the empty-room variance; model vs measured empty room 115 vs 115 fT, 21.4 vs 20.2 fT/cm.
- Conditions: intrinsic, intrinsic+brain, intrinsic+brain+env, projected. Detectability = known-topography matched-filter SNR sqrt(s^T C^+ s); not an event detection rate or a localization accuracy.

## Neuromag magnetometers (102)

- Geometry: 102 channels, 1 mag + 2 planar grad per site; scalp-to-sensor distance median 29.8 mm (5-95 %: 25.4-35.9 mm).
- Retained rank: intrinsic 102, intrinsic+brain 102, intrinsic+brain+env 102, projected 99 (channel subset after the full-array projection).
- Noise composition (median RMS): intrinsic 20.7 fT, brain 192.1 fT, room 113.0 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.37 | 1.29 | 0.98 | 0.67 | 0.46 | 0.33 | 0.26 | 0.22 | 0.17 | 0.15 | 0.12 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.35 | 1.28 | 0.96 | 0.66 | 0.44 | 0.31 | 0.24 | 0.19 | 0.15 | 0.12 | 0.10 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.95 (10 s), 0.99 (60 s), projected: 0.93 (10 s), 0.98 (60 s).

## Neuromag planar gradiometers (204)

- Geometry: 204 channels, 1 mag + 2 planar grad per site; scalp-to-sensor distance median 29.8 mm (5-95 %: 25.4-35.9 mm).
- Retained rank: intrinsic 204, intrinsic+brain 204, intrinsic+brain+env 204, projected 204 (channel subset after the full-array projection).
- Noise composition (median RMS): intrinsic 21.3 fT/cm, brain 37.1 fT/cm, room 1.7 fT/cm.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.30 | 1.21 | 0.92 | 0.63 | 0.44 | 0.30 | 0.24 | 0.20 | 0.15 | 0.13 | 0.11 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.30 | 1.21 | 0.92 | 0.63 | 0.43 | 0.30 | 0.24 | 0.20 | 0.15 | 0.12 | 0.10 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.89 (10 s), 0.98 (60 s), projected: 0.90 (10 s), 0.98 (60 s).

## Neuromag combined (306)

- Geometry: 306 channels, 1 mag + 2 planar grad per site; scalp-to-sensor distance median 29.8 mm (5-95 %: 25.4-35.9 mm).
- Retained rank: intrinsic 306, intrinsic+brain 306, intrinsic+brain+env 306, projected 298.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.48 | 1.38 | 1.03 | 0.69 | 0.47 | 0.34 | 0.27 | 0.23 | 0.17 | 0.15 | 0.13 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.48 | 1.37 | 1.02 | 0.68 | 0.47 | 0.32 | 0.25 | 0.21 | 0.16 | 0.14 | 0.11 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.85 (10 s), 0.97 (60 s), projected: 0.87 (10 s), 0.97 (60 s).

## OPM matched sites

- Geometry: 98 channels, 98 sites, 1 (scalp normal); scalp-to-sensor distance median 6.9 mm (5-95 %: 6.3-7.0 mm); role: matched-site coverage control.
- Retained rank: intrinsic 98, intrinsic+brain 98, intrinsic+brain+env 98, projected 90.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 524.8 fT, room 116.3 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.72 | 1.57 | 1.10 | 0.71 | 0.47 | 0.33 | 0.26 | 0.21 | 0.16 | 0.14 | 0.12 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.69 | 1.53 | 1.05 | 0.65 | 0.40 | 0.25 | 0.17 | 0.13 | 0.09 | 0.08 | 0.07 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.53x [0.51-0.55], OPM higher for 0 % of targets, 0 % of parcels | 2.30x [2.24-2.36], OPM higher for 100 % of targets, 100 % of parcels | 0.55x [0.53-0.57], OPM higher for 0 % of targets, 0 % of parcels |
| intrinsic+brain | 1.02x [1.00-1.04], OPM higher for 60 % of targets, 53 % of parcels | 1.14x [1.12-1.16], OPM higher for 92 % of targets, 90 % of parcels | 1.05x [1.03-1.07], OPM higher for 72 % of targets, 64 % of parcels |
| intrinsic+brain+env | 1.02x [0.99-1.03], OPM higher for 57 % of targets, 49 % of parcels | 1.11x [1.09-1.14], OPM higher for 86 % of targets, 87 % of parcels | 1.07x [1.06-1.09], OPM higher for 85 % of targets, 83 % of parcels |
| projected | 0.92x [0.85-0.95], OPM higher for 32 % of targets, 23 % of parcels | 0.99x [0.91-1.04], OPM higher for 48 % of targets, 37 % of parcels | 0.97x [0.90-1.01], OPM higher for 44 % of targets, 33 % of parcels |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.95 (10 s), 0.99 (60 s), projected: 0.95 (10 s), 0.99 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.02x, 10 mm 1.02x, 20 mm 1.02x.
- Joint OPM noise x scalp gap (vs combined, intrinsic+brain; the arrays are rebuilt for each gap): gap0mm, asd15fT 1.02x; gap0mm, asd20fT 0.99x; gap0mm, asd30fT 0.93x; gap3mm, asd15fT 0.99x; gap3mm, asd20fT 0.96x; gap3mm, asd30fT 0.89x; gap6mm, asd15fT 0.96x; gap6mm, asd20fT 0.92x; gap6mm, asd30fT 0.84x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 8.0 fT/sqrt(Hz), vs grad 34.4 fT/sqrt(Hz), vs mag 8.3 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.10x; opm asd 10fT 1.06x; opm asd 15fT 1.02x; opm asd 20fT 0.99x; opm asd 30fT 0.93x; background corr 5mm 1.02x; background corr 10mm 1.03x; background mag calibrated 1.02x; head x+5mm 1.01x; head x-5mm 1.02x; head y+5mm 1.02x; head y-5mm 1.01x; head z+5mm 1.00x; head z-5mm 1.04x; head pitch+5deg 1.01x; head pitch-5deg 1.02x; head well fitted 1.00x; gap 3mm 0.99x; gap 6mm 0.96x.
- By lobe (vs combined, intrinsic+brain): frontal 1.04x, parietal 1.04x, temporal 0.99x, occipital 1.02x, cingulate 1.00x, insula 0.96x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.02x, projected 0.93x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.02x, 1-10Hz 1.01x, 8-30Hz 1.02x, 30-80Hz 0.99x.

## OPM channel-budget control (204)

- Geometry: 204 channels, 204 sites, 1 (scalp normal); scalp-to-sensor distance median 6.9 mm (5-95 %: 6.2-7.0 mm); role: channel-budget control vs 204 gradiometers.
- Retained rank: intrinsic 204, intrinsic+brain 204, intrinsic+brain+env 204, projected 196.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 500.1 fT, room 115.0 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.40 | 1.92 | 1.28 | 0.79 | 0.53 | 0.38 | 0.30 | 0.25 | 0.21 | 0.19 | 0.15 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.37 | 1.90 | 1.25 | 0.75 | 0.47 | 0.32 | 0.24 | 0.20 | 0.17 | 0.13 | 0.10 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.75x [0.73-0.78], OPM higher for 4 % of targets, 1 % of parcels | 3.24x [3.14-3.32], OPM higher for 100 % of targets, 100 % of parcels | 0.78x [0.75-0.80], OPM higher for 6 % of targets, 3 % of parcels |
| intrinsic+brain | 1.20x [1.18-1.21], OPM higher for 100 % of targets, 100 % of parcels | 1.39x [1.35-1.42], OPM higher for 100 % of targets, 100 % of parcels | 1.25x [1.23-1.27], OPM higher for 100 % of targets, 100 % of parcels |
| intrinsic+brain+env | 1.20x [1.18-1.22], OPM higher for 99 % of targets, 100 % of parcels | 1.36x [1.33-1.39], OPM higher for 100 % of targets, 100 % of parcels | 1.29x [1.27-1.31], OPM higher for 100 % of targets, 100 % of parcels |
| projected | 1.14x [1.11-1.17], OPM higher for 80 % of targets, 89 % of parcels | 1.27x [1.23-1.30], OPM higher for 90 % of targets, 96 % of parcels | 1.22x [1.19-1.24], OPM higher for 89 % of targets, 96 % of parcels |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.91 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.20x, 10 mm 1.20x, 20 mm 1.17x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 11.3 fT/sqrt(Hz), vs grad 48.6 fT/sqrt(Hz), vs mag 11.6 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.35x; opm asd 10fT 1.28x; opm asd 15fT 1.20x; opm asd 20fT 1.13x; opm asd 30fT 1.04x; background corr 5mm 1.21x; background corr 10mm 1.21x; background mag calibrated 1.20x; head x+5mm 1.18x; head x-5mm 1.20x; head y+5mm 1.21x; head y-5mm 1.18x; head z+5mm 1.16x; head z-5mm 1.23x; head pitch+5deg 1.18x; head pitch-5deg 1.21x; head well fitted 1.16x; gap 3mm 1.14x; gap 6mm 1.08x.
- By lobe (vs combined, intrinsic+brain): frontal 1.25x, parietal 1.16x, temporal 1.19x, occipital 1.15x, cingulate 1.20x, insula 1.17x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.20x, projected 1.14x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.20x, 1-10Hz 1.20x, 8-30Hz 1.19x, 30-80Hz 1.14x.

## OPM dense array (full system)

- Geometry: 215 channels, 215 sites, 1 (scalp normal); scalp-to-sensor distance median 6.9 mm (5-95 %: 6.2-7.0 mm); role: densest feasible single-axis array (full system).
- Retained rank: intrinsic 215, intrinsic+brain 215, intrinsic+brain+env 215, projected 207.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 497.7 fT, room 115.4 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.44 | 1.94 | 1.29 | 0.80 | 0.53 | 0.38 | 0.30 | 0.25 | 0.21 | 0.19 | 0.15 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.40 | 1.92 | 1.26 | 0.76 | 0.47 | 0.32 | 0.24 | 0.20 | 0.17 | 0.13 | 0.11 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.77x [0.75-0.80], OPM higher for 5 % of targets, 3 % of parcels | 3.32x [3.23-3.41], OPM higher for 100 % of targets, 100 % of parcels | 0.79x [0.76-0.82], OPM higher for 8 % of targets, 3 % of parcels |
| intrinsic+brain | 1.21x [1.19-1.23], OPM higher for 100 % of targets, 100 % of parcels | 1.40x [1.36-1.44], OPM higher for 100 % of targets, 100 % of parcels | 1.26x [1.24-1.29], OPM higher for 100 % of targets, 100 % of parcels |
| intrinsic+brain+env | 1.21x [1.19-1.23], OPM higher for 99 % of targets, 100 % of parcels | 1.37x [1.34-1.40], OPM higher for 100 % of targets, 100 % of parcels | 1.30x [1.27-1.33], OPM higher for 100 % of targets, 100 % of parcels |
| projected | 1.15x [1.12-1.18], OPM higher for 81 % of targets, 90 % of parcels | 1.28x [1.24-1.31], OPM higher for 91 % of targets, 96 % of parcels | 1.23x [1.20-1.26], OPM higher for 90 % of targets, 97 % of parcels |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.91 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.21x, 10 mm 1.20x, 20 mm 1.18x.
- Joint OPM noise x scalp gap (vs combined, intrinsic+brain; the arrays are rebuilt for each gap): gap0mm, asd15fT 1.21x; gap0mm, asd20fT 1.14x; gap0mm, asd30fT 1.05x; gap3mm, asd15fT 1.15x; gap3mm, asd20fT 1.08x; gap3mm, asd30fT 1.00x; gap6mm, asd15fT 1.09x; gap6mm, asd20fT 1.03x; gap6mm, asd30fT 0.95x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 11.6 fT/sqrt(Hz), vs grad 49.8 fT/sqrt(Hz), vs mag 11.9 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.36x; opm asd 10fT 1.29x; opm asd 15fT 1.21x; opm asd 20fT 1.14x; opm asd 30fT 1.05x; background corr 5mm 1.21x; background corr 10mm 1.22x; background mag calibrated 1.21x; head x+5mm 1.19x; head x-5mm 1.21x; head y+5mm 1.22x; head y-5mm 1.19x; head z+5mm 1.17x; head z-5mm 1.24x; head pitch+5deg 1.19x; head pitch-5deg 1.22x; head well fitted 1.17x; gap 3mm 1.15x; gap 6mm 1.09x.
- By lobe (vs combined, intrinsic+brain): frontal 1.26x, parietal 1.16x, temporal 1.19x, occipital 1.17x, cingulate 1.21x, insula 1.17x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.20x, projected 1.15x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.21x, 1-10Hz 1.21x, 8-30Hz 1.20x, 30-80Hz 1.15x.

## Convergence

- Background grid vs every usable vertex: median log2 ratios change by <= 0.006.
- BEM 5,120 vs 20,480 triangles (1 layer): <= 0.0023; 3 vs 1 layer: <= 0.119 (dense/combined 1.21x with 3 layers, 1.13x with 1 layer).
- Oct-6 vs random full-resolution targets: <= 0.032; Neuromag 4-point vs accurate integration: 0.6 % in the gains.

## Limitations

- Bootstrap CIs resample Desikan-Killiany parcels of one anatomy (targets within a parcel are correlated); they do not include model or between-subject uncertainty. Sensitivity analyses vary one factor at a time except the joint OPM noise x scalp gap grid; they show the dependence, they do not bound it.
- The 'projected' condition removes the simulated room field exactly (it lies in the removed 8-dim subspace); it measures the projection's cost (rank, signal attenuation), not residual interference.
- Head-position variants move the head in the fixed Neuromag helmet; the room field is kept in head coordinates (8-term model; the change over a 5-mm move is second order).
- Detectability is a known-topography matched-filter SNR, not an event detection rate or localization accuracy.
- One adult anatomy and one measured head position; between-subject variability is not represented.
- OPM intrinsic noise is a declared sweep, not a device specification; OPM movement artefacts, cross-talk and calibration errors are not modelled.
- Scalp-gap variants rebuild the OPM arrays (sites farther out pack more easily: 215, 223 and 231 dense sites at 0, 3 and 6 mm), so they mix the gap with extra sensors (about 1 %).
- Known sensor-side BEM error (v2 fix planned): some OPM cell integration points lie within 1 mm of, or inside, the 3-layer BEM head surface; single channels err by up to 72 %, the headline ratios by <= 0.8 % (methods section 8).
