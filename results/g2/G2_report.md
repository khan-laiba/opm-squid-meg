# G2 report: realistic adult OPM vs Neuromag comparison (NEW)

Generated from `results/g2/g2_summary.json` (code commit 19a8fd2; band supplement `g2_band_sensitivity.json` at 19a8fd2; MNE 1.13.2). Methods: `docs/methods.md` section 8; assumptions in `docs/provenance_register.md`. This is a proposed study; no author of the reproduced papers has reviewed it.

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

- Geometry: 95 channels, 95 sites, 1 (scalp normal); scalp-to-sensor distance median 7.7 mm (5-95 %: 6.9-9.0 mm); role: matched-site coverage control.
- Retained rank: intrinsic 95, intrinsic+brain 95, intrinsic+brain+env 95, projected 87.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 508.6 fT, room 116.3 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.70 | 1.54 | 1.09 | 0.70 | 0.46 | 0.32 | 0.25 | 0.20 | 0.15 | 0.14 | 0.11 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.67 | 1.50 | 1.03 | 0.64 | 0.39 | 0.24 | 0.16 | 0.11 | 0.08 | 0.06 | 0.05 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.51x [0.49-0.53], OPM higher for 0 % of targets, 0 % of parcels | 2.21x [2.15-2.27], OPM higher for 100 % of targets, 100 % of parcels | 0.53x [0.51-0.55], OPM higher for 0 % of targets, 0 % of parcels |
| intrinsic+brain | 1.00x [0.97-1.02], OPM higher for 49 % of targets, 39 % of parcels | 1.10x [1.08-1.12], OPM higher for 88 % of targets, 90 % of parcels | 1.02x [1.00-1.05], OPM higher for 60 % of targets, 50 % of parcels |
| intrinsic+brain+env | 0.99x [0.96-1.01], OPM higher for 46 % of targets, 34 % of parcels | 1.08x [1.04-1.11], OPM higher for 75 % of targets, 73 % of parcels | 1.04x [1.02-1.06], OPM higher for 72 % of targets, 67 % of parcels |
| projected | 0.89x [0.80-0.93], OPM higher for 29 % of targets, 17 % of parcels | 0.96x [0.87-1.02], OPM higher for 44 % of targets, 34 % of parcels | 0.93x [0.85-0.98], OPM higher for 39 % of targets, 27 % of parcels |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.95 (10 s), 0.99 (60 s), projected: 0.95 (10 s), 0.99 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.00x, 10 mm 1.00x, 20 mm 1.00x.
- Joint OPM noise x scalp gap (vs combined, intrinsic+brain; the same sites moved outward): gap0mm, asd15fT 1.00x; gap0mm, asd20fT 0.97x; gap0mm, asd30fT 0.92x; gap3mm, asd15fT 0.98x; gap3mm, asd20fT 0.94x; gap3mm, asd30fT 0.87x; gap6mm, asd15fT 0.95x; gap6mm, asd20fT 0.90x; gap6mm, asd30fT 0.82x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 7.7 fT/sqrt(Hz), vs grad 33.1 fT/sqrt(Hz), vs mag 7.9 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.05x; opm asd 10fT 1.03x; opm asd 15fT 1.00x; opm asd 20fT 0.97x; opm asd 30fT 0.92x; background corr 5mm 1.00x; background corr 10mm 1.01x; background mag calibrated 0.99x; head x+5mm 0.99x; head x-5mm 1.00x; head y+5mm 1.00x; head y-5mm 0.99x; head z+5mm 0.98x; head z-5mm 1.01x; head pitch+5deg 0.99x; head pitch-5deg 1.00x; head well fitted 0.98x; gap 3mm 0.98x; gap 6mm 0.95x.
- By lobe (vs combined, intrinsic+brain): frontal 1.03x, parietal 1.03x, temporal 0.94x, occipital 1.00x, cingulate 0.97x, insula 0.91x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.00x, projected 0.91x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.00x, 1-10Hz 0.99x, 8-30Hz 1.00x, 30-80Hz 0.98x.

## OPM channel-budget control (204)

- Geometry: 204 channels, 204 sites, 1 (scalp normal); scalp-to-sensor distance median 7.8 mm (5-95 %: 6.8-9.0 mm); role: channel-budget control vs 204 gradiometers.
- Retained rank: intrinsic 204, intrinsic+brain 204, intrinsic+brain+env 204, projected 196.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 491.8 fT, room 115.0 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.38 | 1.89 | 1.27 | 0.78 | 0.50 | 0.35 | 0.27 | 0.23 | 0.17 | 0.15 | 0.12 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.35 | 1.86 | 1.22 | 0.73 | 0.45 | 0.29 | 0.21 | 0.16 | 0.11 | 0.09 | 0.08 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.74x [0.71-0.76], OPM higher for 3 % of targets, 0 % of parcels | 3.16x [3.08-3.25], OPM higher for 100 % of targets, 100 % of parcels | 0.76x [0.73-0.79], OPM higher for 4 % of targets, 1 % of parcels |
| intrinsic+brain | 1.11x [1.08-1.14], OPM higher for 90 % of targets, 94 % of parcels | 1.28x [1.25-1.31], OPM higher for 100 % of targets, 100 % of parcels | 1.16x [1.11-1.19], OPM higher for 95 % of targets, 99 % of parcels |
| intrinsic+brain+env | 1.11x [1.08-1.15], OPM higher for 85 % of targets, 90 % of parcels | 1.26x [1.21-1.29], OPM higher for 98 % of targets, 100 % of parcels | 1.19x [1.15-1.22], OPM higher for 99 % of targets, 100 % of parcels |
| projected | 1.05x [0.98-1.10], OPM higher for 57 % of targets, 51 % of parcels | 1.15x [1.07-1.21], OPM higher for 69 % of targets, 66 % of parcels | 1.10x [1.04-1.16], OPM higher for 66 % of targets, 61 % of parcels |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.91 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.11x, 10 mm 1.12x, 20 mm 1.11x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 11.0 fT/sqrt(Hz), vs grad 47.4 fT/sqrt(Hz), vs mag 11.4 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.24x; opm asd 10fT 1.18x; opm asd 15fT 1.11x; opm asd 20fT 1.06x; opm asd 30fT 1.00x; background corr 5mm 1.12x; background corr 10mm 1.14x; background mag calibrated 1.11x; head x+5mm 1.10x; head x-5mm 1.11x; head y+5mm 1.12x; head y-5mm 1.10x; head z+5mm 1.08x; head z-5mm 1.14x; head pitch+5deg 1.09x; head pitch-5deg 1.12x; head well fitted 1.08x; gap 3mm 1.07x; gap 6mm 1.02x.
- By lobe (vs combined, intrinsic+brain): frontal 1.21x, parietal 1.13x, temporal 1.10x, occipital 1.07x, cingulate 1.01x, insula 1.04x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.12x, projected 1.07x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.11x, 1-10Hz 1.11x, 8-30Hz 1.11x, 30-80Hz 1.09x.

## OPM dense array (full system)

- Geometry: 205 channels, 205 sites, 1 (scalp normal); scalp-to-sensor distance median 7.8 mm (5-95 %: 6.8-9.0 mm); role: densest feasible single-axis array (full system).
- Retained rank: intrinsic 205, intrinsic+brain 205, intrinsic+brain+env 205, projected 197.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 490.9 fT, room 114.7 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.38 | 1.89 | 1.27 | 0.78 | 0.50 | 0.35 | 0.27 | 0.23 | 0.17 | 0.15 | 0.12 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.35 | 1.86 | 1.22 | 0.73 | 0.45 | 0.29 | 0.21 | 0.16 | 0.11 | 0.09 | 0.08 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.74x [0.71-0.76], OPM higher for 3 % of targets, 0 % of parcels | 3.17x [3.09-3.25], OPM higher for 100 % of targets, 100 % of parcels | 0.76x [0.73-0.79], OPM higher for 4 % of targets, 1 % of parcels |
| intrinsic+brain | 1.11x [1.08-1.14], OPM higher for 90 % of targets, 94 % of parcels | 1.28x [1.25-1.32], OPM higher for 100 % of targets, 100 % of parcels | 1.16x [1.11-1.19], OPM higher for 95 % of targets, 99 % of parcels |
| intrinsic+brain+env | 1.11x [1.08-1.15], OPM higher for 85 % of targets, 90 % of parcels | 1.26x [1.21-1.29], OPM higher for 98 % of targets, 100 % of parcels | 1.19x [1.16-1.22], OPM higher for 99 % of targets, 100 % of parcels |
| projected | 1.05x [0.99-1.10], OPM higher for 57 % of targets, 51 % of parcels | 1.15x [1.07-1.21], OPM higher for 69 % of targets, 66 % of parcels | 1.10x [1.04-1.16], OPM higher for 66 % of targets, 61 % of parcels |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.91 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.12x, 10 mm 1.12x, 20 mm 1.11x.
- Joint OPM noise x scalp gap (vs combined, intrinsic+brain; the same sites moved outward): gap0mm, asd15fT 1.11x; gap0mm, asd20fT 1.06x; gap0mm, asd30fT 1.00x; gap3mm, asd15fT 1.07x; gap3mm, asd20fT 1.02x; gap3mm, asd30fT 0.96x; gap6mm, asd15fT 1.02x; gap6mm, asd20fT 0.98x; gap6mm, asd30fT 0.91x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 11.1 fT/sqrt(Hz), vs grad 47.5 fT/sqrt(Hz), vs mag 11.4 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.24x; opm asd 10fT 1.18x; opm asd 15fT 1.11x; opm asd 20fT 1.06x; opm asd 30fT 1.00x; background corr 5mm 1.12x; background corr 10mm 1.14x; background mag calibrated 1.11x; head x+5mm 1.10x; head x-5mm 1.11x; head y+5mm 1.12x; head y-5mm 1.10x; head z+5mm 1.08x; head z-5mm 1.14x; head pitch+5deg 1.09x; head pitch-5deg 1.12x; head well fitted 1.08x; gap 3mm 1.07x; gap 6mm 1.02x.
- By lobe (vs combined, intrinsic+brain): frontal 1.21x, parietal 1.13x, temporal 1.10x, occipital 1.07x, cingulate 1.01x, insula 1.04x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.12x, projected 1.07x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.11x, 1-10Hz 1.11x, 8-30Hz 1.11x, 30-80Hz 1.09x.

## Convergence

- Background grid vs every usable vertex: median log2 ratios change by <= 0.007.
- 3-layer head surface 20,480 triangles (primary, A-BEM-SKIN) vs the v1 5,120: dense/combined 1.109x vs 1.110x, matched/combined 0.996x vs 0.996x (intrinsic+brain, convergence subset).
- BEM 5,120 vs 20,480 triangles (1 layer): <= 0.0017; 3 vs 1 layer: <= 0.027 (dense/combined 1.11x with 3 layers, 1.11x with 1 layer).
- Oct-6 vs random full-resolution targets: <= 0.038; Neuromag 4-point vs accurate integration: 0.6 % in the gains.

## Limitations

- Bootstrap CIs resample Desikan-Killiany parcels of one anatomy (targets within a parcel are correlated); they do not include model or between-subject uncertainty. Sensitivity analyses vary one factor at a time except the joint OPM noise x scalp gap grid; they show the dependence, they do not bound it.
- The 'projected' condition removes the simulated room field exactly (it lies in the removed 8-dim subspace); it measures the projection's cost (rank, signal attenuation), not residual interference.
- Head-position variants move the head in the fixed Neuromag helmet; the room field is kept in head coordinates (8-term model; the change over a 5-mm move is second order).
- Detectability is a known-topography matched-filter SNR, not an event detection rate or localization accuracy.
- One adult anatomy and one measured head position; between-subject variability is not represented.
- OPM intrinsic noise is a declared sweep, not a device specification; OPM movement artefacts, cross-talk and calibration errors are not modelled.
- Head model: 3-layer BEM with the head surface refined to 20,480 triangles and the whole OPM cell >= 1 mm outside it at the sampled points (v3; exact cube-to-mesh distance >= 1.001 mm). Refining the head surface changes the headline by -0.1 %; if the error falls with the square of the mesh size (not verified on this head), the refined surface is within ~0.03 % (methods section 3). The 1-layer model gives nearly the same dense/combined ratio (convergence section).
- Scalp-gap variants move the primary OPM sites outward along their axes (same sites).
