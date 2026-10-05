# G2 report: realistic adult OPM vs Neuromag comparison (NEW)

Generated from `results/g2/g2_summary.json` (code commit e53bea8; band supplement `g2_band_sensitivity.json` at e53bea8; MNE 1.13.2). Methods: `docs/methods.md` section 8; assumptions in `docs/provenance_register.md`. This is a proposed study; no author of the reproduced papers has reviewed it.

## Common setup

- Anatomy: MNE sample subject, measured head position; 7661 target dipoles (10 nAm, cortical normal, usable oct-6 vertices); background grid of 1755 area-weighted sources.
- Neuromag geometry: the sample recording's Vectorview sensor positions and transforms with MRN T3 coil types (3024 magnetometers, 3014 planar gradiometers), a representative Neuromag system, not one installation; the manufacturer's typical noise values (MEGIN TRIUX; not verified against the specification sheet).
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
| median detectability (10 nAm) | 1.35 | 1.28 | 0.96 | 0.66 | 0.44 | 0.31 | 0.24 | 0.19 | 0.15 | 0.11 | 0.10 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.95 (10 s), 0.99 (60 s), projected: 0.93 (10 s), 0.98 (60 s).

## Neuromag planar gradiometers (204)

- Geometry: 204 channels, 1 mag + 2 planar grad per site; scalp-to-sensor distance median 29.8 mm (5-95 %: 25.4-35.9 mm).
- Retained rank: intrinsic 204, intrinsic+brain 204, intrinsic+brain+env 204, projected 204 (channel subset after the full-array projection).
- Noise composition (median RMS): intrinsic 21.3 fT/cm, brain 37.1 fT/cm, room 1.7 fT/cm.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.30 | 1.21 | 0.92 | 0.63 | 0.44 | 0.30 | 0.24 | 0.20 | 0.15 | 0.13 | 0.10 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.30 | 1.21 | 0.92 | 0.63 | 0.43 | 0.30 | 0.24 | 0.19 | 0.15 | 0.12 | 0.10 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.89 (10 s), 0.98 (60 s), projected: 0.89 (10 s), 0.98 (60 s).

## Neuromag combined (306)

- Geometry: 306 channels, 1 mag + 2 planar grad per site; scalp-to-sensor distance median 29.8 mm (5-95 %: 25.4-35.9 mm).
- Retained rank: intrinsic 306, intrinsic+brain 306, intrinsic+brain+env 306, projected 298.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.48 | 1.38 | 1.03 | 0.69 | 0.47 | 0.34 | 0.27 | 0.22 | 0.17 | 0.15 | 0.13 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.48 | 1.38 | 1.02 | 0.69 | 0.46 | 0.32 | 0.25 | 0.21 | 0.16 | 0.13 | 0.11 |

- Estimated (plug-in) covariance, detectability relative to oracle: intrinsic+brain: 0.85 (10 s), 0.97 (60 s), projected: 0.86 (10 s), 0.97 (60 s).

## OPM matched sites

- Geometry: 98 channels, 98 sites, 1 (scalp normal); scalp-to-sensor distance median 7.0 mm (5-95 %: 6.9-7.1 mm); role: matched-site coverage control.
- Retained rank: intrinsic 98, intrinsic+brain 98, intrinsic+brain+env 98, projected 90.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 524.9 fT, room 115.8 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.72 | 1.57 | 1.10 | 0.71 | 0.47 | 0.33 | 0.26 | 0.21 | 0.16 | 0.14 | 0.12 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 1.70 | 1.53 | 1.06 | 0.67 | 0.43 | 0.29 | 0.21 | 0.16 | 0.12 | 0.09 | 0.08 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.53x [0.51-0.55], OPM higher for 0 % of targets, 0 % of parcels | 2.29x [2.23-2.34], OPM higher for 100 % of targets, 100 % of parcels | 0.55x [0.53-0.57], OPM higher for 0 % of targets, 0 % of parcels |
| intrinsic+brain | 1.01x [0.99-1.03], OPM higher for 56 % of targets, 47 % of parcels | 1.13x [1.11-1.15], OPM higher for 92 % of targets, 90 % of parcels | 1.04x [1.01-1.06], OPM higher for 70 % of targets, 61 % of parcels |
| intrinsic+brain+env | 1.01x [0.99-1.03], OPM higher for 54 % of targets, 47 % of parcels | 1.11x [1.08-1.13], OPM higher for 86 % of targets, 87 % of parcels | 1.06x [1.05-1.08], OPM higher for 86 % of targets, 81 % of parcels |
| projected | 0.95x [0.90-0.98], OPM higher for 36 % of targets, 27 % of parcels | 1.02x [0.97-1.07], OPM higher for 55 % of targets, 50 % of parcels | 1.00x [0.96-1.03], OPM higher for 49 % of targets, 41 % of parcels |

- Estimated (plug-in) covariance (from one noise realization shared by all arrays), detectability relative to oracle: intrinsic+brain: 0.95 (10 s), 0.99 (60 s), projected: 0.95 (10 s), 0.99 (60 s).
- Metric dependence (vs combined, intrinsic+brain): detectability 1.01x; peak-channel SNR (best single channel) 0.85x [0.83-0.87], the OPM array higher for 17 % of targets (Neuromag ahead); mean-power SNR 1.05x.
- After the external-field projection, by depth (vs combined): 10 mm 1.19x, 15 mm 1.11x, 20 mm 1.04x, 25 mm 0.99x, 30 mm 0.93x, 35 mm 0.88x, 40 mm 0.83x, 45 mm 0.82x, 50 mm 0.79x, 55 mm 0.75x, 60 mm 0.76x; the OPM array is behind Neuromag from 25 mm down.
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.01x, 10 mm 1.01x, 20 mm 1.01x.
- Joint OPM noise x scalp gap (vs combined, intrinsic+brain; the same sites moved outward): gap0mm, asd15fT 1.01x; gap0mm, asd20fT 0.98x; gap0mm, asd30fT 0.93x; gap3mm, asd15fT 0.99x; gap3mm, asd20fT 0.96x; gap3mm, asd30fT 0.89x; gap6mm, asd15fT 0.96x; gap6mm, asd20fT 0.91x; gap6mm, asd30fT 0.84x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 8.0 fT/sqrt(Hz), vs grad 34.3 fT/sqrt(Hz), vs mag 8.2 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.08x; opm asd 10fT 1.05x; opm asd 15fT 1.01x; opm asd 20fT 0.98x; opm asd 30fT 0.93x; background corr 5mm 1.01x; background corr 10mm 1.02x; background mag calibrated 1.01x; squid measured spectrum 1.02x; head x+5mm 1.00x; head x-5mm 1.01x; head y+5mm 1.01x; head y-5mm 1.00x; head z+5mm 0.99x; head z-5mm 1.02x; head pitch+5deg 1.00x; head pitch-5deg 1.01x; head well fitted 0.99x; gap 3mm 0.99x; gap 6mm 0.96x.
- By lobe (vs combined, intrinsic+brain): frontal 1.04x, parietal 1.04x, temporal 0.99x, occipital 1.02x, cingulate 0.99x, insula 0.95x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.02x, projected 0.96x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.01x, 1-10Hz 1.00x, 8-30Hz 1.01x, 30-80Hz 0.99x.

## OPM channel-budget control (204)

- Geometry: 204 channels, 204 sites, 1 (scalp normal); scalp-to-sensor distance median 7.0 mm (5-95 %: 6.8-7.5 mm); role: channel-budget control vs 204 gradiometers.
- Retained rank: intrinsic 204, intrinsic+brain 204, intrinsic+brain+env 204, projected 196.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 494.4 fT, room 114.2 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.43 | 1.90 | 1.27 | 0.79 | 0.51 | 0.37 | 0.28 | 0.24 | 0.18 | 0.16 | 0.13 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.41 | 1.89 | 1.25 | 0.76 | 0.49 | 0.33 | 0.25 | 0.20 | 0.15 | 0.13 | 0.10 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.73x [0.71-0.76], OPM higher for 3 % of targets, 1 % of parcels | 3.14x [3.06-3.22], OPM higher for 100 % of targets, 100 % of parcels | 0.76x [0.73-0.78], OPM higher for 5 % of targets, 3 % of parcels |
| intrinsic+brain | 1.14x [1.12-1.16], OPM higher for 100 % of targets, 100 % of parcels | 1.33x [1.29-1.35], OPM higher for 100 % of targets, 100 % of parcels | 1.19x [1.16-1.22], OPM higher for 100 % of targets, 100 % of parcels |
| intrinsic+brain+env | 1.14x [1.12-1.16], OPM higher for 99 % of targets, 100 % of parcels | 1.29x [1.26-1.32], OPM higher for 100 % of targets, 100 % of parcels | 1.22x [1.20-1.25], OPM higher for 100 % of targets, 100 % of parcels |
| projected | 1.11x [1.08-1.14], OPM higher for 83 % of targets, 87 % of parcels | 1.24x [1.21-1.28], OPM higher for 93 % of targets, 97 % of parcels | 1.20x [1.17-1.23], OPM higher for 94 % of targets, 96 % of parcels |

- Estimated (plug-in) covariance (from one noise realization shared by all arrays), detectability relative to oracle: intrinsic+brain: 0.91 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Metric dependence (vs combined, intrinsic+brain): detectability 1.14x; peak-channel SNR (best single channel) 0.90x [0.88-0.92], the OPM array higher for 28 % of targets (Neuromag ahead); mean-power SNR 1.04x.
- After the external-field projection, by depth (vs combined): 10 mm 1.60x, 15 mm 1.36x, 20 mm 1.21x, 25 mm 1.11x, 30 mm 1.05x, 35 mm 1.02x, 40 mm 1.01x, 45 mm 1.00x, 50 mm 0.99x, 55 mm 0.96x, 60 mm 0.92x; the OPM array is behind Neuromag from 45 mm down.
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.14x, 10 mm 1.15x, 20 mm 1.13x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 11.0 fT/sqrt(Hz), vs grad 47.1 fT/sqrt(Hz), vs mag 11.4 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.29x; opm asd 10fT 1.22x; opm asd 15fT 1.14x; opm asd 20fT 1.08x; opm asd 30fT 1.01x; background corr 5mm 1.15x; background corr 10mm 1.17x; background mag calibrated 1.14x; squid measured spectrum 1.16x; head x+5mm 1.13x; head x-5mm 1.14x; head y+5mm 1.15x; head y-5mm 1.12x; head z+5mm 1.11x; head z-5mm 1.17x; head pitch+5deg 1.13x; head pitch-5deg 1.15x; head well fitted 1.10x; gap 3mm 1.09x; gap 6mm 1.04x.
- By lobe (vs combined, intrinsic+brain): frontal 1.21x, parietal 1.14x, temporal 1.15x, occipital 1.11x, cingulate 1.05x, insula 1.08x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.14x, projected 1.12x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.14x, 1-10Hz 1.14x, 8-30Hz 1.14x, 30-80Hz 1.10x.

## OPM dense array (full system)

- Geometry: 208 channels, 208 sites, 1 (scalp normal); scalp-to-sensor distance median 7.0 mm (5-95 %: 6.8-7.5 mm); role: densest feasible single-axis array (full system); the densest array found under the 17-mm centre-spacing rule, not proven maximal; the package footprint is an unverified assumption (U-OPM-PACK).
- Retained rank: intrinsic 208, intrinsic+brain 208, intrinsic+brain+env 208, projected 200.
- Noise composition (median RMS): intrinsic 88.8 fT, brain 494.4 fT, room 114.3 fT.

Detectability vs depth, intrinsic+brain:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.43 | 1.91 | 1.28 | 0.79 | 0.51 | 0.37 | 0.28 | 0.24 | 0.18 | 0.16 | 0.13 |

Detectability vs depth, projected:

| depth [mm] | 10-15 | 15-20 | 20-25 | 25-30 | 30-35 | 35-40 | 40-45 | 45-50 | 50-55 | 55-60 | 60-65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median detectability (10 nAm) | 2.41 | 1.90 | 1.25 | 0.76 | 0.49 | 0.33 | 0.25 | 0.20 | 0.15 | 0.13 | 0.10 |

OPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):

| condition | vs combined | vs gradiometers | vs magnetometers |
|---|---|---|---|
| intrinsic | 0.74x [0.72-0.77], OPM higher for 4 % of targets, 1 % of parcels | 3.17x [3.08-3.25], OPM higher for 100 % of targets, 100 % of parcels | 0.76x [0.74-0.79], OPM higher for 6 % of targets, 3 % of parcels |
| intrinsic+brain | 1.14x [1.12-1.17], OPM higher for 100 % of targets, 100 % of parcels | 1.33x [1.30-1.36], OPM higher for 100 % of targets, 100 % of parcels | 1.19x [1.17-1.22], OPM higher for 100 % of targets, 100 % of parcels |
| intrinsic+brain+env | 1.14x [1.12-1.17], OPM higher for 99 % of targets, 100 % of parcels | 1.29x [1.26-1.32], OPM higher for 100 % of targets, 100 % of parcels | 1.23x [1.20-1.26], OPM higher for 100 % of targets, 100 % of parcels |
| projected | 1.12x [1.08-1.15], OPM higher for 84 % of targets, 87 % of parcels | 1.25x [1.21-1.29], OPM higher for 94 % of targets, 97 % of parcels | 1.20x [1.17-1.23], OPM higher for 95 % of targets, 96 % of parcels |

- Estimated (plug-in) covariance (from one noise realization shared by all arrays), detectability relative to oracle: intrinsic+brain: 0.91 (10 s), 0.98 (60 s), projected: 0.91 (10 s), 0.98 (60 s).
- Metric dependence (vs combined, intrinsic+brain): detectability 1.14x; peak-channel SNR (best single channel) 0.90x [0.88-0.93], the OPM array higher for 28 % of targets (Neuromag ahead); mean-power SNR 1.04x.
- After the external-field projection, by depth (vs combined): 10 mm 1.61x, 15 mm 1.36x, 20 mm 1.21x, 25 mm 1.11x, 30 mm 1.05x, 35 mm 1.03x, 40 mm 1.01x, 45 mm 1.00x, 50 mm 0.99x, 55 mm 0.96x, 60 mm 0.92x; the OPM array is behind Neuromag from 45 mm down.
- Patches vs Neuromag combined (intrinsic+brain): 5 mm 1.14x, 10 mm 1.15x, 20 mm 1.13x.
- Joint OPM noise x scalp gap (vs combined, intrinsic+brain; the same sites moved outward): gap0mm, asd15fT 1.14x; gap0mm, asd20fT 1.08x; gap0mm, asd30fT 1.01x; gap3mm, asd15fT 1.09x; gap3mm, asd20fT 1.03x; gap3mm, asd30fT 0.97x; gap6mm, asd15fT 1.04x; gap6mm, asd20fT 0.99x; gap6mm, asd30fT 0.92x.
- Intrinsic noise only, break-even OPM noise (ratio = 1): vs combined 11.1 fT/sqrt(Hz), vs grad 47.6 fT/sqrt(Hz), vs mag 11.5 fT/sqrt(Hz).
- Sensitivity, one factor at a time (vs combined, intrinsic+brain): opm asd 7fT 1.30x; opm asd 10fT 1.23x; opm asd 15fT 1.14x; opm asd 20fT 1.08x; opm asd 30fT 1.01x; background corr 5mm 1.15x; background corr 10mm 1.17x; background mag calibrated 1.15x; squid measured spectrum 1.17x; head x+5mm 1.13x; head x-5mm 1.14x; head y+5mm 1.16x; head y-5mm 1.13x; head z+5mm 1.11x; head z-5mm 1.18x; head pitch+5deg 1.13x; head pitch-5deg 1.15x; head well fitted 1.10x; gap 3mm 1.09x; gap 6mm 1.04x.
- By lobe (vs combined, intrinsic+brain): frontal 1.22x, parietal 1.14x, temporal 1.16x, occipital 1.12x, cingulate 1.05x, insula 1.09x.
- Without the 558 medial-wall targets (7.3 %; FreeSurfer 'unknown', not cortex), vs combined: intrinsic+brain 1.15x, projected 1.13x.
- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): 1-40Hz 1.14x, 1-10Hz 1.15x, 8-30Hz 1.14x, 30-80Hz 1.11x.

## Channel-count control: a triaxial OPM at the matched sites (A-OPM-TRIAX)

98 matched sites x 3 axes = 294 channels (Neuromag: 102 sites x 3 = 306): the normal and the two tangential axes of each 10-mm cell, 15 fT/sqrt(Hz) on every axis, or the tangential axes at twice that. Median detectability ratio (95 % CI, parcel bootstrap):

| tangential noise | condition | vs combined | vs gradiometers | vs magnetometers | triaxial / normal axis only |
|---|---|---|---|---|---|
| equal | intrinsic+brain | 1.12x [1.10-1.15], OPM higher for 95 % of targets, 94 % of parcels | 1.30x [1.27-1.33], OPM higher for 100 % of targets, 100 % of parcels | 1.17x [1.14-1.20], OPM higher for 100 % of targets, 100 % of parcels | 1.12x [1.10-1.13], OPM higher for 100 % of targets, 100 % of parcels |
| equal | projected | 1.17x [1.15-1.19], OPM higher for 99 % of targets, 100 % of parcels | 1.30x [1.27-1.34], OPM higher for 100 % of targets, 100 % of parcels | 1.26x [1.23-1.28], OPM higher for 100 % of targets, 100 % of parcels | 1.23x [1.20-1.28], OPM higher for 100 % of targets, 100 % of parcels |
| x2 | intrinsic+brain | 1.05x [1.03-1.08], OPM higher for 78 % of targets, 70 % of parcels | 1.21x [1.18-1.23], OPM higher for 99 % of targets, 100 % of parcels | 1.09x [1.06-1.12], OPM higher for 90 % of targets, 87 % of parcels | 1.05x [1.04-1.06], OPM higher for 100 % of targets, 100 % of parcels |
| x2 | projected | 1.09x [1.08-1.11], OPM higher for 93 % of targets, 89 % of parcels | 1.21x [1.18-1.24], OPM higher for 99 % of targets, 100 % of parcels | 1.17x [1.15-1.19], OPM higher for 98 % of targets, 100 % of parcels | 1.15x [1.13-1.19], OPM higher for 100 % of targets, 100 % of parcels |

## Neuromag sensor noise from the measured spectrum

- Per channel, the in-band empty-room variance the 8-term room fit leaves (an upper bound on the sensor noise in this room): median 25.9 fT (magnetometers) and 20.0 fT/cm (gradiometers), vs the brochure 20.7 fT and 21.3 fT/cm in the band. Dense / matched vs combined: intrinsic+brain 1.17x / 1.02x, projected 1.13x / 0.96x.

## Link to the analytical benchmark (G1A)

- The sphere model with the real standoffs (OPM 7.0 mm, Neuromag 29.8 mm median; peak-field ratio, eta = 3) gives an equal-SNR depth of 29.0 mm (Jas's 0 / 18 mm standoffs: 27.7 mm); the matched OPM array on this head, by its peak-field ratio, 29.6 mm; it is ahead at every depth for eta <= 2.25 and behind at every depth for eta >= 5. The realistic comparison above replaces eta by explicit noise and the peak field by the known-topography detectability of all channels (methods section 8).

## Convergence

- Background grid vs every usable vertex: median log2 ratios change by <= 0.005.
- 3-layer head surface 20,480 triangles (primary, A-BEM-SKIN) vs the v1 5,120: dense/combined 1.141x vs 1.149x, matched/combined 1.008x vs 1.008x (intrinsic+brain, convergence subset).
- BEM 5,120 vs 20,480 triangles (1 layer): <= 0.0015; 3 vs 1 layer: <= 0.046 (dense/combined 1.14x with 3 layers, 1.13x with 1 layer).
- Oct-6 vs random full-resolution targets: <= 0.029; Neuromag 4-point vs accurate integration: 0.6 % in the gains.

## Limitations

- Bootstrap CIs resample Desikan-Killiany parcels of one anatomy (targets within a parcel are correlated); they do not include model or between-subject uncertainty. Sensitivity analyses vary one factor at a time except the joint OPM noise x scalp gap grid; they show the dependence, they do not bound it.
- The 'projected' condition removes the simulated room field exactly (it lies in the removed 8-dim subspace); it measures the projection's cost (rank, signal attenuation), not residual interference.
- Head-position variants move the head in the fixed Neuromag helmet; the room field is kept in head coordinates (8-term model; the change over a 5-mm move is second order).
- Detectability is a known-topography matched-filter SNR, not an event detection rate or localization accuracy.
- One adult anatomy and one measured head position; between-subject variability is not represented.
- OPM intrinsic noise is a declared sweep, not a device specification; OPM movement artefacts, cross-talk and calibration errors are not modelled.
- Head model: 3-layer BEM, the head surface on the MRI scalp (A-BEM-CONFORM) and refined to 20,480 triangles; the whole OPM cell >= 1 mm outside it at the sampled points (exact cube-to-mesh distance: matched 1.056 mm, dense 1.006 mm). Refining the head surface from 5,120 to 20,480 triangles changes the dense/combined headline from 1.149x to 1.141x on a 1,000-target subset (`bem_skin_refinement.json`, computed at ed852b2; this report's G2 results: e53bea8). The 1-layer model gives nearly the same dense/combined ratio (convergence section).
- Scalp-gap variants move the primary OPM sites outward along their axes (same sites).
