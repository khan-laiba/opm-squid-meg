# G2 report: realistic adult OPM vs Neuromag comparison (NEW)

Generated from `results/g2/g2_summary.json` (code commit 26985d4; band supplement `g2_band_sensitivity.json` at 26985d4; MNE 1.13.2). Methods: `docs/methods.md` section 8; assumptions in `docs/provenance_register.md`. This is a proposed study; no author of the reproduced papers has reviewed it.

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
| median detectability (10 nAm) | 1.48 | 1.37 | 1.02 | 0.68 | 0.46 | 0.32 | 0.25 | 0.21 | 0.16 | 0.14 | 0.11 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.85 (10 s), 0.97 (60 s), projected: 0.87 (10 s), 0.97 (60 s).

## OPM matched sites

- Geometry: 97 channels, 97 sites, 1 (scalp normal); scalp-to-sensor distance median 7.0 mm (5-95 %: 6.4-7.8 mm); role: matched-site coverage control.
- Retained rank: intrinsic 97, intrinsic+brain 97, intrinsic+brain+env 97, projected 89.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 526.1 fT, room 116.3 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.72 | 1.57 | 1.10 | 0.71 | 0.47 | 0.32 | 0.25 | 0.20 | 0.15 | 0.14 | 0.11 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.69 | 1.53 | 1.04 | 0.64 | 0.39 | 0.25 | 0.17 | 0.12 | 0.08 | 0.07 | 0.06 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.53x [0.51-0.55], OPM higher for 0 % of targets, 0 % of parcels | 2.28x [2.21-2.34], OPM higher for 100 % of targets, 100 % of parcels | 0.55x [0.52-0.57], OPM higher for 0 % of targets, 0 % of parcels |
| intrinsic+brain | 1.00x [0.98-1.03], OPM higher for 52 % of targets, 44 % of parcels | 1.12x [1.10-1.14], OPM higher for 90 % of targets, 89 % of parcels | 1.03x [1.00-1.06], OPM higher for 64 % of targets, 54 % of parcels |
| intrinsic+brain+env | 1.00x [0.97-1.02], OPM higher for 50 % of targets, 41 % of parcels | 1.09x [1.06-1.12], OPM higher for 80 % of targets, 76 % of parcels | 1.05x [1.04-1.07], OPM higher for 78 % of targets, 73 % of parcels |
| projected | 0.90x [0.81-0.95], OPM higher for 31 % of targets, 20 % of parcels | 0.97x [0.88-1.03], OPM higher for 46 % of targets, 37 % of parcels | 0.95x [0.87-1.00], OPM higher for 41 % of targets, 33 % of parcels |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.95 (10 s), 0.99 (60 s), projected: 0.95 (10 s), 0.99 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.01x, 10 mm 1.01x, 20 mm 1.01x.
- Joint OPM noise x scalp gap (vs combined, intrinsic+brain; the arrays are rebuilt for each gap): gap0mm, asd15fT 1.00x; gap0mm, asd20fT 0.98x; gap0mm, asd30fT 0.93x; gap3mm, asd15fT 0.98x; gap3mm, asd20fT 0.95x; gap3mm, asd30fT 0.88x; gap6mm, asd15fT 0.96x; gap6mm, asd20fT 0.91x; gap6mm, asd30fT 0.83x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 7.9 fT/sqrt(Hz), vs grad 34.2 fT/sqrt(Hz), vs mag 8.2 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.07x; opm asd 10fT 1.04x; opm asd 15fT 1.00x; opm asd 20fT 0.98x; opm asd 30fT 0.93x; background corr 5mm 1.01x; background corr 10mm 1.02x; background mag calibrated 1.00x; head x+5mm 1.00x; head x-5mm 1.01x; head y+5mm 1.01x; head y-5mm 1.00x; head z+5mm 0.99x; head z-5mm 1.02x; head pitch+5deg 1.00x; head pitch-5deg 1.01x; head well fitted 0.98x; gap 3mm 0.98x; gap 6mm 0.96x.
- By lobe (vs combined, intrinsic+brain): frontal 1.04x, parietal 1.04x, temporal 0.96x, occipital 1.01x, cingulate 0.98x, insula 0.93x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.01x, projected 0.92x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.00x, 1-10Hz 1.00x, 8-30Hz 1.01x, 30-80Hz 0.99x.

## OPM channel-budget control (204)

- Geometry: 204 channels, 204 sites, 1 (scalp normal); scalp-to-sensor distance median 6.9 mm (5-95 %: 6.5-7.8 mm); role: channel-budget control vs 204 gradiometers.
- Retained rank: intrinsic 204, intrinsic+brain 204, intrinsic+brain+env 204, projected 196.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 502.5 fT, room 114.6 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.41 | 1.93 | 1.28 | 0.78 | 0.51 | 0.36 | 0.28 | 0.23 | 0.17 | 0.15 | 0.13 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.38 | 1.90 | 1.24 | 0.74 | 0.46 | 0.30 | 0.22 | 0.16 | 0.12 | 0.10 | 0.08 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.75x [0.73-0.78], OPM higher for 4 % of targets, 1 % of parcels | 3.24x [3.15-3.33], OPM higher for 100 % of targets, 100 % of parcels | 0.78x [0.75-0.80], OPM higher for 7 % of targets, 3 % of parcels |
| intrinsic+brain | 1.12x [1.09-1.16], OPM higher for 96 % of targets, 100 % of parcels | 1.30x [1.27-1.34], OPM higher for 100 % of targets, 100 % of parcels | 1.17x [1.14-1.21], OPM higher for 99 % of targets, 100 % of parcels |
| intrinsic+brain+env | 1.13x [1.09-1.16], OPM higher for 91 % of targets, 96 % of parcels | 1.27x [1.24-1.31], OPM higher for 99 % of targets, 100 % of parcels | 1.21x [1.18-1.24], OPM higher for 100 % of targets, 100 % of parcels |
| projected | 1.07x [1.01-1.12], OPM higher for 61 % of targets, 54 % of parcels | 1.17x [1.10-1.23], OPM higher for 73 % of targets, 70 % of parcels | 1.13x [1.07-1.18], OPM higher for 71 % of targets, 69 % of parcels |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.91 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.13x, 10 mm 1.14x, 20 mm 1.12x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 11.3 fT/sqrt(Hz), vs grad 48.6 fT/sqrt(Hz), vs mag 11.6 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.26x; opm asd 10fT 1.20x; opm asd 15fT 1.12x; opm asd 20fT 1.07x; opm asd 30fT 1.00x; background corr 5mm 1.13x; background corr 10mm 1.15x; background mag calibrated 1.12x; head x+5mm 1.11x; head x-5mm 1.12x; head y+5mm 1.13x; head y-5mm 1.11x; head z+5mm 1.09x; head z-5mm 1.15x; head pitch+5deg 1.11x; head pitch-5deg 1.13x; head well fitted 1.09x; gap 3mm 1.08x; gap 6mm 1.03x.
- By lobe (vs combined, intrinsic+brain): frontal 1.22x, parietal 1.14x, temporal 1.13x, occipital 1.07x, cingulate 1.02x, insula 1.06x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.14x, projected 1.09x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.12x, 1-10Hz 1.12x, 8-30Hz 1.12x, 30-80Hz 1.10x.

## OPM dense array (full system)

- Geometry: 211 channels, 211 sites, 1 (scalp normal); scalp-to-sensor distance median 6.9 mm (5-95 %: 6.5-7.7 mm); role: densest feasible single-axis array (full system).
- Retained rank: intrinsic 211, intrinsic+brain 211, intrinsic+brain+env 211, projected 203.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 499.3 fT, room 114.6 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.43 | 1.94 | 1.29 | 0.79 | 0.51 | 0.36 | 0.28 | 0.23 | 0.17 | 0.15 | 0.13 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.40 | 1.91 | 1.25 | 0.75 | 0.46 | 0.30 | 0.22 | 0.16 | 0.12 | 0.10 | 0.08 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.76x [0.74-0.79], OPM higher for 5 % of targets, 1 % of parcels | 3.29x [3.20-3.37], OPM higher for 100 % of targets, 100 % of parcels | 0.79x [0.76-0.82], OPM higher for 8 % of targets, 3 % of parcels |
| intrinsic+brain | 1.13x [1.10-1.16], OPM higher for 96 % of targets, 100 % of parcels | 1.31x [1.28-1.35], OPM higher for 100 % of targets, 100 % of parcels | 1.18x [1.14-1.21], OPM higher for 99 % of targets, 100 % of parcels |
| intrinsic+brain+env | 1.13x [1.10-1.17], OPM higher for 91 % of targets, 97 % of parcels | 1.28x [1.24-1.32], OPM higher for 99 % of targets, 100 % of parcels | 1.21x [1.18-1.25], OPM higher for 100 % of targets, 100 % of parcels |
| projected | 1.07x [1.01-1.12], OPM higher for 62 % of targets, 56 % of parcels | 1.18x [1.10-1.24], OPM higher for 73 % of targets, 71 % of parcels | 1.14x [1.08-1.18], OPM higher for 71 % of targets, 69 % of parcels |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.91 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.13x, 10 mm 1.14x, 20 mm 1.13x.
- Joint OPM noise x scalp gap (vs combined, intrinsic+brain; the arrays are rebuilt for each gap): gap0mm, asd15fT 1.13x; gap0mm, asd20fT 1.08x; gap0mm, asd30fT 1.01x; gap3mm, asd15fT 1.08x; gap3mm, asd20fT 1.03x; gap3mm, asd30fT 0.97x; gap6mm, asd15fT 1.03x; gap6mm, asd20fT 0.99x; gap6mm, asd30fT 0.93x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 11.5 fT/sqrt(Hz), vs grad 49.4 fT/sqrt(Hz), vs mag 11.8 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.27x; opm asd 10fT 1.20x; opm asd 15fT 1.13x; opm asd 20fT 1.08x; opm asd 30fT 1.01x; background corr 5mm 1.14x; background corr 10mm 1.16x; background mag calibrated 1.13x; head x+5mm 1.11x; head x-5mm 1.13x; head y+5mm 1.14x; head y-5mm 1.11x; head z+5mm 1.10x; head z-5mm 1.16x; head pitch+5deg 1.11x; head pitch-5deg 1.14x; head well fitted 1.09x; gap 3mm 1.08x; gap 6mm 1.03x.
- By lobe (vs combined, intrinsic+brain): frontal 1.23x, parietal 1.14x, temporal 1.13x, occipital 1.08x, cingulate 1.03x, insula 1.06x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.14x, projected 1.09x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.13x, 1-10Hz 1.13x, 8-30Hz 1.13x, 30-80Hz 1.11x.

## Convergence

- Background grid vs every usable vertex: median log2 ratios change by <= 0.007.
- 3-layer head surface 20,480 triangles (primary, A-BEM-SKIN) vs the v1 5,120: dense/combined 1.127x vs 1.132x, matched/combined 1.004x vs 1.007x (intrinsic+brain, convergence subset).
- BEM 5,120 vs 20,480 triangles (1 layer): <= 0.0043; 3 vs 1 layer: <= 0.029 (dense/combined 1.13x with 3 layers, 1.12x with 1 layer).
- Oct-6 vs random full-resolution targets: <= 0.036; Neuromag 4-point vs accurate integration: 0.6 % in the gains.

## Limitations

- Bootstrap CIs resample Desikan-Killiany parcels of one anatomy (targets within a parcel are correlated); they do not include model or between-subject uncertainty. Sensitivity analyses vary one factor at a time except the joint OPM noise x scalp gap grid; they show the dependence, they do not bound it.
- The 'projected' condition removes the simulated room field exactly (it lies in the removed 8-dim subspace); it measures the projection's cost (rank, signal attenuation), not residual interference.
- Head-position variants move the head in the fixed Neuromag helmet; the room field is kept in head coordinates (8-term model; the change over a 5-mm move is second order).
- Detectability is a known-topography matched-filter SNR, not an event detection rate or localization accuracy.
- One adult anatomy and one measured head position; between-subject variability is not represented.
- OPM intrinsic noise is a declared sweep, not a device specification; OPM movement artefacts, cross-talk and calibration errors are not modelled.
- Head model: 3-layer BEM with the head surface refined to 20,480 triangles and every OPM cell integration point >= 1 mm outside it (v2). Near the head surface the BEM field is only approximately converged (about 1 % for the headline, methods section 3); the 1-layer model gives a lower dense/combined ratio (convergence section).
- Scalp-gap variants move the primary OPM sites outward along their axes (same sites).
