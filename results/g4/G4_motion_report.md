# G4 extension: head motion and OPM slippage (NEW, bounded secondary analysis)

Static fit was studied first (G2, G3B, G4). This analysis adds what the static maps cannot show, on the G3B arrays and noise conventions (Neuromag at top contact; the refitted dense OPM array; intrinsic + brain noise, the G3B primary condition; 10-nAm cortical-normal dipoles; area-weighted medians over cortical targets). It bounds two mechanisms; it does not establish motion robustness. Configuration: `configs/g4_motion.toml`; code: `scripts/g4_motion.py`, `src/opmsquid/motion.py`.

## A. Sustained displacement: the head in the fixed helmet (Neuromag) vs a slipped cap (dense OPM)

'known': the displaced geometry is known (template and noise of the displaced geometry: ideal movement compensation or a known slip). 'mismatched': the template of the reference geometry applied to the displaced data (noise covariance from those data). A head-mounted array moving with the head has no geometry change; only slips appear for it. A wrong-polarity template counts as -60 dB.

| anatomy | system | displacement | known [dB] | mismatched [dB] | share > 3 dB loss (mismatched) | note |
|---|---|---|---|---|---|---|
| adult | Neuromag combined | down 2 mm | -0.08 | -0.13 | 0.00 | nearest magnetometer 22.0 mm |
| adult | Neuromag combined | x+2 mm | +0.01 | -0.05 | 0.00 | nearest magnetometer 19.5 mm |
| adult | Neuromag combined | x-2 mm | -0.00 | -0.06 | 0.00 | nearest magnetometer 20.1 mm |
| adult | Neuromag combined | y+2 mm | -0.01 | -0.07 | 0.00 | nearest magnetometer 20.6 mm |
| adult | Neuromag combined | y-2 mm | +0.02 | -0.04 | 0.00 | nearest magnetometer 19.7 mm |
| adult | Neuromag combined | down 5 mm | -0.19 | -0.53 | 0.00 | nearest magnetometer 23.0 mm |
| adult | Neuromag combined | x+5 mm | +0.02 | -0.33 | 0.00 | nearest magnetometer 18.2 mm |
| adult | Neuromag combined | x-5 mm | -0.00 | -0.38 | 0.00 | nearest magnetometer 19.0 mm |
| adult | Neuromag combined | y+5 mm | -0.03 | -0.35 | 0.00 | nearest magnetometer 18.9 mm |
| adult | Neuromag combined | y-5 mm | +0.04 | -0.28 | 0.00 | nearest magnetometer 18.7 mm |
| adult | Neuromag combined | down 10 mm | -0.38 | -1.65 | 0.05 | nearest magnetometer 23.2 mm |
| adult | Neuromag | x+10 mm | infeasible | | | nearest magnetometer 14.0 mm |
| adult | Neuromag | x-10 mm | infeasible | | | nearest magnetometer 14.5 mm |
| adult | Neuromag | y+10 mm | infeasible | | | nearest magnetometer 16.6 mm |
| adult | Neuromag | y-10 mm | infeasible | | | nearest magnetometer 16.2 mm |
| adult | Neuromag combined | pitch +5 deg | +0.02 | -0.68 | 0.04 | nearest magnetometer 19.2 mm |
| adult | Neuromag combined | pitch -5 deg | -0.02 | -0.84 | 0.04 | nearest magnetometer 19.8 mm |
| adult | Neuromag combined | roll +5 deg | +0.00 | -0.69 | 0.06 | nearest magnetometer 19.0 mm |
| adult | Neuromag | roll -5 deg | infeasible | | | nearest magnetometer 17.6 mm |
| adult | Neuromag combined | yaw +5 deg | -0.00 | -0.30 | 0.00 | nearest magnetometer 20.2 mm |
| adult | Neuromag combined | yaw -5 deg | +0.00 | -0.29 | 0.00 | nearest magnetometer 20.3 mm |
| adult | Neuromag | pitch +10 deg | infeasible | | | nearest magnetometer 16.8 mm |
| adult | Neuromag | pitch -10 deg | infeasible | | | nearest magnetometer 15.6 mm |
| adult | Neuromag | roll +10 deg | infeasible | | | nearest magnetometer 16.3 mm |
| adult | Neuromag | roll -10 deg | infeasible | | | nearest magnetometer 12.5 mm |
| adult | Neuromag combined | yaw +10 deg | -0.00 | -1.20 | 0.15 | nearest magnetometer 19.9 mm |
| adult | Neuromag combined | yaw -10 deg | +0.01 | -1.19 | 0.14 | nearest magnetometer 18.2 mm |
| adult | OPM dense | slip x +1 deg | -0.03 | -0.09 | 0.00 | sensors moved 1.9 mm (median); 35 lifted (max 3.0 mm) |
| adult | OPM dense | slip x -1 deg | +0.01 | -0.05 | 0.00 | sensors moved 1.9 mm (median); 26 lifted (max 1.5 mm) |
| adult | OPM dense | slip y +1 deg | -0.02 | -0.09 | 0.00 | sensors moved 1.7 mm (median); 53 lifted (max 1.5 mm) |
| adult | OPM dense | slip y -1 deg | -0.02 | -0.07 | 0.00 | sensors moved 1.7 mm (median); 50 lifted (max 1.0 mm) |
| adult | OPM dense | slip z +1 deg | -0.01 | -0.05 | 0.00 | sensors moved 1.5 mm (median); 27 lifted (max 1.0 mm) |
| adult | OPM dense | slip z -1 deg | -0.00 | -0.04 | 0.00 | sensors moved 1.5 mm (median); 22 lifted (max 1.5 mm) |
| adult | OPM dense | slip x +3 deg | -0.08 | -0.53 | 0.02 | sensors moved 5.7 mm (median); 77 lifted (max 9.0 mm) |
| adult | OPM dense | slip x -3 deg | -0.02 | -0.43 | 0.01 | sensors moved 5.5 mm (median); 86 lifted (max 24.0 mm) |
| adult | OPM dense | slip y +3 deg | -0.08 | -0.47 | 0.01 | sensors moved 5.0 mm (median); 87 lifted (max 3.5 mm) |
| adult | OPM dense | slip y -3 deg | -0.05 | -0.43 | 0.01 | sensors moved 5.0 mm (median); 79 lifted (max 3.5 mm) |
| adult | OPM dense | slip z +3 deg | -0.02 | -0.23 | 0.00 | sensors moved 4.5 mm (median); 55 lifted (max 3.5 mm) |
| adult | OPM dense | slip z -3 deg | -0.01 | -0.21 | 0.00 | sensors moved 4.5 mm (median); 46 lifted (max 3.5 mm) |
| 2-year template | Neuromag combined | down 2 mm | -0.07 | -0.13 | 0.00 | nearest magnetometer 21.6 mm |
| 2-year template | Neuromag combined | x+2 mm | +0.02 | -0.04 | 0.00 | nearest magnetometer 18.9 mm |
| 2-year template | Neuromag combined | x-2 mm | -0.02 | -0.08 | 0.00 | nearest magnetometer 21.1 mm |
| 2-year template | Neuromag combined | y+2 mm | -0.01 | -0.07 | 0.00 | nearest magnetometer 20.8 mm |
| 2-year template | Neuromag combined | y-2 mm | +0.01 | -0.04 | 0.00 | nearest magnetometer 19.6 mm |
| 2-year template | Neuromag combined | down 5 mm | -0.18 | -0.54 | 0.00 | nearest magnetometer 23.4 mm |
| 2-year template | Neuromag | x+5 mm | infeasible | | | nearest magnetometer 16.7 mm |
| 2-year template | Neuromag combined | x-5 mm | -0.03 | -0.43 | 0.00 | nearest magnetometer 22.5 mm |
| 2-year template | Neuromag combined | y+5 mm | -0.02 | -0.37 | 0.00 | nearest magnetometer 21.1 mm |
| 2-year template | Neuromag combined | y-5 mm | +0.04 | -0.29 | 0.00 | nearest magnetometer 18.9 mm |
| 2-year template | Neuromag combined | down 10 mm | -0.36 | -1.67 | 0.07 | nearest magnetometer 25.5 mm |
| 2-year template | Neuromag | x+10 mm | infeasible | | | nearest magnetometer 12.9 mm |
| 2-year template | Neuromag combined | x-10 mm | -0.04 | -1.46 | 0.07 | nearest magnetometer 22.6 mm |
| 2-year template | Neuromag combined | y+10 mm | -0.03 | -1.37 | 0.06 | nearest magnetometer 21.0 mm |
| 2-year template | Neuromag | y-10 mm | infeasible | | | nearest magnetometer 17.3 mm |
| 2-year template | Neuromag combined | pitch +5 deg | +0.03 | -0.48 | 0.02 | nearest magnetometer 19.7 mm |
| 2-year template | Neuromag combined | pitch -5 deg | -0.03 | -0.62 | 0.02 | nearest magnetometer 20.4 mm |
| 2-year template | Neuromag combined | roll +5 deg | -0.00 | -0.50 | 0.04 | nearest magnetometer 21.5 mm |
| 2-year template | Neuromag | roll -5 deg | infeasible | | | nearest magnetometer 17.6 mm |
| 2-year template | Neuromag combined | yaw +5 deg | -0.00 | -0.26 | 0.00 | nearest magnetometer 20.1 mm |
| 2-year template | Neuromag combined | yaw -5 deg | +0.00 | -0.24 | 0.00 | nearest magnetometer 20.3 mm |
| 2-year template | Neuromag combined | pitch +10 deg | +0.06 | -2.04 | 0.38 | nearest magnetometer 19.3 mm |
| 2-year template | Neuromag combined | pitch -10 deg | -0.05 | -2.31 | 0.41 | nearest magnetometer 19.7 mm |
| 2-year template | Neuromag combined | roll +10 deg | +0.01 | -1.95 | 0.35 | nearest magnetometer 22.4 mm |
| 2-year template | Neuromag | roll -10 deg | infeasible | | | nearest magnetometer 14.9 mm |
| 2-year template | Neuromag combined | yaw +10 deg | -0.00 | -1.04 | 0.10 | nearest magnetometer 19.7 mm |
| 2-year template | Neuromag combined | yaw -10 deg | +0.01 | -0.98 | 0.10 | nearest magnetometer 20.6 mm |
| 2-year template | OPM dense | slip x +1 deg | -0.01 | -0.05 | 0.00 | sensors moved 1.6 mm (median); 1 lifted (max 0.5 mm) |
| 2-year template | OPM dense | slip x -1 deg | +0.01 | -0.02 | 0.00 | sensors moved 1.6 mm (median); 0 lifted (max 0.0 mm) |
| 2-year template | OPM dense | slip y +1 deg | -0.00 | -0.03 | 0.00 | sensors moved 1.5 mm (median); 5 lifted (max 0.5 mm) |
| 2-year template | OPM dense | slip y -1 deg | +0.00 | -0.03 | 0.00 | sensors moved 1.5 mm (median); 1 lifted (max 0.5 mm) |
| 2-year template | OPM dense | slip z +1 deg | +0.00 | -0.01 | 0.00 | sensors moved 1.3 mm (median); 0 lifted (max 0.0 mm) |
| 2-year template | OPM dense | slip z -1 deg | -0.00 | -0.01 | 0.00 | sensors moved 1.3 mm (median); 1 lifted (max 0.5 mm) |
| 2-year template | OPM dense | slip x +3 deg | -0.05 | -0.38 | 0.02 | sensors moved 4.9 mm (median); 18 lifted (max 3.5 mm) |
| 2-year template | OPM dense | slip x -3 deg | -0.00 | -0.28 | 0.01 | sensors moved 4.7 mm (median); 49 lifted (max 2.0 mm) |
| 2-year template | OPM dense | slip y +3 deg | -0.03 | -0.31 | 0.02 | sensors moved 4.2 mm (median); 43 lifted (max 2.5 mm) |
| 2-year template | OPM dense | slip y -3 deg | -0.03 | -0.31 | 0.02 | sensors moved 4.3 mm (median); 44 lifted (max 2.5 mm) |
| 2-year template | OPM dense | slip z +3 deg | +0.00 | -0.11 | 0.00 | sensors moved 3.9 mm (median); 6 lifted (max 1.0 mm) |
| 2-year template | OPM dense | slip z -3 deg | -0.00 | -0.12 | 0.00 | sensors moved 3.9 mm (median); 10 lifted (max 1.5 mm) |
| 12-month template | Neuromag combined | down 2 mm | -0.08 | -0.14 | 0.00 | nearest magnetometer 21.8 mm |
| 12-month template | Neuromag combined | x+2 mm | +0.01 | -0.05 | 0.00 | nearest magnetometer 19.2 mm |
| 12-month template | Neuromag combined | x-2 mm | -0.01 | -0.07 | 0.00 | nearest magnetometer 20.9 mm |
| 12-month template | Neuromag combined | y+2 mm | -0.01 | -0.07 | 0.00 | nearest magnetometer 20.9 mm |
| 12-month template | Neuromag combined | y-2 mm | +0.01 | -0.05 | 0.00 | nearest magnetometer 19.6 mm |
| 12-month template | Neuromag combined | down 5 mm | -0.20 | -0.56 | 0.00 | nearest magnetometer 24.3 mm |
| 12-month template | Neuromag | x+5 mm | infeasible | | | nearest magnetometer 17.8 mm |
| 12-month template | Neuromag combined | x-5 mm | -0.02 | -0.41 | 0.00 | nearest magnetometer 20.7 mm |
| 12-month template | Neuromag combined | y+5 mm | -0.02 | -0.38 | 0.00 | nearest magnetometer 20.7 mm |
| 12-month template | Neuromag combined | y-5 mm | +0.03 | -0.33 | 0.00 | nearest magnetometer 18.7 mm |
| 12-month template | Neuromag combined | down 10 mm | -0.39 | -1.72 | 0.08 | nearest magnetometer 28.4 mm |
| 12-month template | Neuromag | x+10 mm | infeasible | | | nearest magnetometer 15.7 mm |
| 12-month template | Neuromag combined | x-10 mm | -0.01 | -1.43 | 0.07 | nearest magnetometer 20.0 mm |
| 12-month template | Neuromag combined | y+10 mm | -0.02 | -1.41 | 0.06 | nearest magnetometer 20.6 mm |
| 12-month template | Neuromag | y-10 mm | infeasible | | | nearest magnetometer 16.8 mm |
| 12-month template | Neuromag combined | pitch +5 deg | +0.03 | -0.47 | 0.02 | nearest magnetometer 19.8 mm |
| 12-month template | Neuromag combined | pitch -5 deg | -0.03 | -0.59 | 0.02 | nearest magnetometer 19.8 mm |
| 12-month template | Neuromag combined | roll +5 deg | +0.00 | -0.49 | 0.04 | nearest magnetometer 20.7 mm |
| 12-month template | Neuromag combined | roll -5 deg | +0.00 | -0.48 | 0.03 | nearest magnetometer 19.0 mm |
| 12-month template | Neuromag combined | yaw +5 deg | -0.00 | -0.24 | 0.00 | nearest magnetometer 20.1 mm |
| 12-month template | Neuromag combined | yaw -5 deg | +0.00 | -0.22 | 0.00 | nearest magnetometer 20.4 mm |
| 12-month template | Neuromag combined | pitch +10 deg | +0.06 | -1.96 | 0.36 | nearest magnetometer 19.5 mm |
| 12-month template | Neuromag combined | pitch -10 deg | -0.06 | -2.25 | 0.40 | nearest magnetometer 19.3 mm |
| 12-month template | Neuromag combined | roll +10 deg | +0.01 | -1.90 | 0.34 | nearest magnetometer 20.5 mm |
| 12-month template | Neuromag | roll -10 deg | infeasible | | | nearest magnetometer 17.2 mm |
| 12-month template | Neuromag combined | yaw +10 deg | -0.00 | -0.94 | 0.06 | nearest magnetometer 19.7 mm |
| 12-month template | Neuromag combined | yaw -10 deg | +0.00 | -0.91 | 0.06 | nearest magnetometer 20.6 mm |
| 12-month template | OPM dense | slip x +1 deg | -0.02 | -0.05 | 0.00 | sensors moved 1.6 mm (median); 3 lifted (max 1.5 mm) |
| 12-month template | OPM dense | slip x -1 deg | +0.01 | -0.02 | 0.00 | sensors moved 1.6 mm (median); 0 lifted (max 0.0 mm) |
| 12-month template | OPM dense | slip y +1 deg | -0.00 | -0.04 | 0.00 | sensors moved 1.4 mm (median); 3 lifted (max 1.5 mm) |
| 12-month template | OPM dense | slip y -1 deg | -0.00 | -0.04 | 0.00 | sensors moved 1.4 mm (median); 3 lifted (max 1.0 mm) |
| 12-month template | OPM dense | slip z +1 deg | +0.00 | -0.01 | 0.00 | sensors moved 1.3 mm (median); 1 lifted (max 0.5 mm) |
| 12-month template | OPM dense | slip z -1 deg | -0.00 | -0.01 | 0.00 | sensors moved 1.3 mm (median); 1 lifted (max 0.5 mm) |
| 12-month template | OPM dense | slip x +3 deg | -0.07 | -0.41 | 0.01 | sensors moved 4.8 mm (median); 27 lifted (max 4.0 mm) |
| 12-month template | OPM dense | slip x -3 deg | +0.00 | -0.31 | 0.02 | sensors moved 4.6 mm (median); 42 lifted (max 2.0 mm) |
| 12-month template | OPM dense | slip y +3 deg | -0.04 | -0.34 | 0.02 | sensors moved 4.2 mm (median); 41 lifted (max 4.0 mm) |
| 12-month template | OPM dense | slip y -3 deg | -0.04 | -0.34 | 0.01 | sensors moved 4.3 mm (median); 41 lifted (max 2.5 mm) |
| 12-month template | OPM dense | slip z +3 deg | -0.00 | -0.13 | 0.00 | sensors moved 3.8 mm (median); 7 lifted (max 1.0 mm) |
| 12-month template | OPM dense | slip z -3 deg | -0.00 | -0.13 | 0.00 | sensors moved 3.8 mm (median); 6 lifted (max 2.0 mm) |

## B. In-band motion of the head-mounted array in a static residual field

The dense OPM array moves rigidly with the head (rotation about a pivot 60 mm below the head origin) in a static field: 1 nT uniform, or a 1 nT/m symmetric traceless gradient (Frobenius norm), random isotropic draws. The artefact covariance is J Sigma J^T for in-band rotation with the given RMS per axis; it is added to the noise after the correction ('none'; 'homogeneous': 3-term homogeneous-field projection; 'homogeneous+gradient': the 8-term projection of G2), applied to signal and noise alike. Calibration errors: RMS tilt of each sensitive axis and RMS gain error, unknown to the analyst. Everything is linear in rotation x field: for a field of B nT the rotation thresholds below divide by B. Neuromag is fixed in the room and has no such term; its static detectability is the comparator for D.

Static reference (no motion): OPM detectability change from the correction alone and D (dense OPM vs Neuromag combined, both intrinsic + brain):

| anatomy | correction | OPM change [dB] | D [dB] |
|---|---|---|---|
| adult | none | +0.00 | +0.99 |
| adult | homogeneous | -0.13 | +0.80 |
| adult | homogeneous+gradient | -0.44 | +0.51 |
| 2-year template | none | +0.00 | +1.85 |
| 2-year template | homogeneous | -0.08 | +1.57 |
| 2-year template | homogeneous+gradient | -0.39 | +1.25 |
| 12-month template | none | +0.00 | +2.06 |
| 12-month template | homogeneous | -0.06 | +1.86 |
| 12-month template | homogeneous+gradient | -0.35 | +1.54 |

In-band artefact per channel (median over channels, then over 32 draws) for 1 deg RMS rotation per axis in the unit field, after each correction [fT]. Then the in-band rotation (deg RMS per axis, in the unit field) at which the median OPM detectability falls by 1 or 3 dB, or D falls to 0 dB, when the artefact is not part of the analyst's noise model (matched filter of the static covariance applied to data that contain it), from the median curve over draws, with the 10th-90th percentiles of the per-draw thresholds ('-': not reached up to 5 deg). The draws are common random numbers: the same fields and calibration errors in every correction, pivot and anatomy. Last column: the loss at 5 deg when the artefact is part of the known noise covariance (the optimal filter nulls its at most 3 spatial patterns: the bound for data-driven nulling or regression on measured head motion).

| anatomy | field | pivot | correction | calibration (tilt, gain) | artefact [fT per deg] | 1 dB | 3 dB | D = 0 | oracle loss at 5 deg [dB] |
|---|---|---|---|---|---|---|---|---|---|
| adult | uniform | neck | none | tilt0deg_gain0pct | 1.5e+04 | 0.0219 [0.0207-0.0243] | 0.0444 [0.0385-0.0532] | 0.0216 [0.0205-0.0236] | -0.06 |
| adult | gradient | neck | none | tilt0deg_gain0pct | 1.83e+03 | 0.169 [0.154-0.193] | 0.316 [0.295-0.347] | 0.155 [0.145-0.172] | -0.09 |
| adult | uniform | neck | none | tilt1deg_gain1pct | 1.5e+04 | 0.0219 [0.0208-0.0242] | 0.0442 [0.0388-0.053] | 0.0215 [0.0206-0.0236] | -0.06 |
| adult | gradient | neck | none | tilt1deg_gain1pct | 1.82e+03 | 0.168 [0.153-0.192] | 0.315 [0.293-0.346] | 0.154 [0.145-0.172] | -0.09 |
| adult | uniform | neck | none | tilt3deg_gain3pct | 1.5e+04 | 0.0217 [0.0206-0.0237] | 0.0433 [0.0381-0.0517] | 0.0214 [0.0205-0.0231] | -0.05 |
| adult | gradient | neck | none | tilt3deg_gain3pct | 1.81e+03 | 0.165 [0.15-0.184] | 0.31 [0.289-0.336] | 0.153 [0.143-0.17] | -0.08 |
| adult | uniform | neck | homogeneous | tilt0deg_gain0pct | 6.07e-12 | - | - | - | +0.00 |
| adult | gradient | neck | homogeneous | tilt0deg_gain0pct | 1.39e+03 | 0.146 [0.133-0.161] | 0.283 [0.262-0.305] | 0.122 [0.115-0.134] | -0.15 |
| adult | uniform | neck | homogeneous | tilt1deg_gain1pct | 236 | 0.449 [0.403-0.507] | 0.898 [0.829-1] | 0.329 [0.303-0.347] | -0.01 |
| adult | gradient | neck | homogeneous | tilt1deg_gain1pct | 1.39e+03 | 0.146 [0.132-0.162] | 0.283 [0.261-0.306] | 0.122 [0.115-0.133] | -0.14 |
| adult | uniform | neck | homogeneous | tilt3deg_gain3pct | 704 | 0.142 [0.129-0.154] | 0.277 [0.255-0.295] | 0.116 [0.11-0.125] | -0.01 |
| adult | gradient | neck | homogeneous | tilt3deg_gain3pct | 1.39e+03 | 0.143 [0.132-0.155] | 0.279 [0.26-0.297] | 0.121 [0.115-0.132] | -0.12 |
| adult | gradient | origin | homogeneous | tilt0deg_gain0pct | 1.39e+03 | 0.146 [0.133-0.161] | 0.283 [0.262-0.305] | 0.122 [0.115-0.134] | -0.15 |
| adult | gradient | origin | homogeneous | tilt1deg_gain1pct | 1.39e+03 | 0.146 [0.132-0.162] | 0.283 [0.261-0.306] | 0.122 [0.115-0.133] | -0.14 |
| adult | gradient | origin | homogeneous | tilt3deg_gain3pct | 1.39e+03 | 0.144 [0.132-0.157] | 0.28 [0.261-0.299] | 0.122 [0.115-0.132] | -0.13 |
| adult | uniform | neck | homogeneous+gradient | tilt0deg_gain0pct | 6.29e-12 | - | - | - | +0.00 |
| adult | gradient | neck | homogeneous+gradient | tilt0deg_gain0pct | 1.27e-12 | - | - | - | +0.00 |
| adult | uniform | neck | homogeneous+gradient | tilt1deg_gain1pct | 234 | 0.382 [0.347-0.453] | 0.795 [0.736-0.905] | 0.232 [0.222-0.248] | -0.01 |
| adult | gradient | neck | homogeneous+gradient | tilt1deg_gain1pct | 27.5 | 2.99 [2.71-3.37] | - | 1.99 [1.73-2.05] | -0.02 |
| adult | uniform | neck | homogeneous+gradient | tilt3deg_gain3pct | 693 | 0.13 [0.118-0.137] | 0.256 [0.235-0.269] | 0.0809 [0.0753-0.0885] | -0.01 |
| adult | gradient | neck | homogeneous+gradient | tilt3deg_gain3pct | 82.9 | 1.11 [1.02-1.19] | 2.21 [2.01-2.38] | 0.613 [0.582-0.658] | -0.02 |
| adult | gradient | origin | homogeneous+gradient | tilt0deg_gain0pct | 1.22e-12 | - | - | - | +0.00 |
| adult | gradient | origin | homogeneous+gradient | tilt1deg_gain1pct | 23.3 | 3.53 [3.05-3.92] | - | 2.18 [2.07-2.27] | -0.02 |
| adult | gradient | origin | homogeneous+gradient | tilt3deg_gain3pct | 71.6 | 1.23 [1.14-1.34] | 2.45 [2.28-2.64] | 0.712 [0.677-0.787] | -0.02 |
| 2-year template | uniform | neck | none | tilt0deg_gain0pct | 1.49e+04 | 0.0213 [0.0178-0.0227] | 0.0411 [0.0328-0.0487] | 0.0264 [0.0229-0.0296] | -0.05 |
| 2-year template | gradient | neck | none | tilt0deg_gain0pct | 1.56e+03 | 0.202 [0.174-0.211] | 0.364 [0.323-0.403] | 0.243 [0.227-0.257] | -0.09 |
| 2-year template | uniform | neck | none | tilt1deg_gain1pct | 1.5e+04 | 0.0213 [0.0178-0.0228] | 0.0412 [0.0328-0.0491] | 0.0264 [0.023-0.0297] | -0.05 |
| 2-year template | gradient | neck | none | tilt1deg_gain1pct | 1.55e+03 | 0.202 [0.173-0.211] | 0.362 [0.322-0.401] | 0.243 [0.227-0.256] | -0.09 |
| 2-year template | uniform | neck | none | tilt3deg_gain3pct | 1.49e+04 | 0.0211 [0.0178-0.0225] | 0.0405 [0.0329-0.0473] | 0.0262 [0.0231-0.0292] | -0.05 |
| 2-year template | gradient | neck | none | tilt3deg_gain3pct | 1.56e+03 | 0.194 [0.173-0.209] | 0.349 [0.322-0.392] | 0.241 [0.225-0.254] | -0.09 |
| 2-year template | uniform | neck | homogeneous | tilt0deg_gain0pct | 6.19e-12 | - | - | - | +0.00 |
| 2-year template | gradient | neck | homogeneous | tilt0deg_gain0pct | 1.2e+03 | 0.199 [0.183-0.208] | 0.354 [0.334-0.392] | 0.231 [0.223-0.244] | -0.14 |
| 2-year template | uniform | neck | homogeneous | tilt1deg_gain1pct | 236 | 0.442 [0.394-0.505] | 0.887 [0.814-0.995] | 0.545 [0.522-0.571] | -0.02 |
| 2-year template | gradient | neck | homogeneous | tilt1deg_gain1pct | 1.2e+03 | 0.199 [0.182-0.209] | 0.354 [0.333-0.393] | 0.232 [0.222-0.245] | -0.14 |
| 2-year template | uniform | neck | homogeneous | tilt3deg_gain3pct | 711 | 0.146 [0.133-0.16] | 0.283 [0.261-0.304] | 0.192 [0.166-0.203] | -0.02 |
| 2-year template | gradient | neck | homogeneous | tilt3deg_gain3pct | 1.19e+03 | 0.194 [0.177-0.206] | 0.349 [0.327-0.381] | 0.23 [0.221-0.243] | -0.13 |
| 2-year template | gradient | origin | homogeneous | tilt0deg_gain0pct | 1.2e+03 | 0.199 [0.183-0.208] | 0.354 [0.334-0.392] | 0.231 [0.223-0.244] | -0.14 |
| 2-year template | gradient | origin | homogeneous | tilt1deg_gain1pct | 1.19e+03 | 0.198 [0.182-0.209] | 0.354 [0.333-0.393] | 0.232 [0.222-0.244] | -0.14 |
| 2-year template | gradient | origin | homogeneous | tilt3deg_gain3pct | 1.18e+03 | 0.197 [0.179-0.207] | 0.353 [0.33-0.383] | 0.231 [0.222-0.243] | -0.13 |
| 2-year template | uniform | neck | homogeneous+gradient | tilt0deg_gain0pct | 6.58e-12 | - | - | - | +0.00 |
| 2-year template | gradient | neck | homogeneous+gradient | tilt0deg_gain0pct | 8.99e-13 | - | - | - | +0.00 |
| 2-year template | uniform | neck | homogeneous+gradient | tilt1deg_gain1pct | 232 | 0.397 [0.368-0.477] | 0.818 [0.772-0.939] | 0.425 [0.383-0.471] | -0.02 |
| 2-year template | gradient | neck | homogeneous+gradient | tilt1deg_gain1pct | 23.3 | 3.92 [3.28-4.49] | - | 3.77 [3.38-4.14] | -0.03 |
| 2-year template | uniform | neck | homogeneous+gradient | tilt3deg_gain3pct | 691 | 0.134 [0.123-0.147] | 0.264 [0.245-0.285] | 0.143 [0.127-0.151] | -0.02 |
| 2-year template | gradient | neck | homogeneous+gradient | tilt3deg_gain3pct | 72.9 | 1.27 [1.2-1.44] | 2.53 [2.39-2.8] | 1.25 [1.18-1.35] | -0.03 |
| 2-year template | gradient | origin | homogeneous+gradient | tilt0deg_gain0pct | 7.66e-13 | - | - | - | +0.00 |
| 2-year template | gradient | origin | homogeneous+gradient | tilt1deg_gain1pct | 20.2 | - | - | 4.81 [4.15-4.84] | -0.03 |
| 2-year template | gradient | origin | homogeneous+gradient | tilt3deg_gain3pct | 61.1 | 1.53 [1.4-1.67] | 2.93 [2.73-3.13] | 1.54 [1.39-1.61] | -0.03 |
| 12-month template | uniform | neck | none | tilt0deg_gain0pct | 1.5e+04 | 0.0228 [0.0204-0.0238] | 0.0489 [0.037-0.0521] | 0.0329 [0.0266-0.0359] | -0.03 |
| 12-month template | gradient | neck | none | tilt0deg_gain0pct | 1.55e+03 | 0.209 [0.189-0.223] | 0.396 [0.342-0.461] | 0.282 [0.26-0.313] | -0.07 |
| 12-month template | uniform | neck | none | tilt1deg_gain1pct | 1.5e+04 | 0.0227 [0.0203-0.0238] | 0.0485 [0.037-0.052] | 0.0328 [0.0266-0.0359] | -0.03 |
| 12-month template | gradient | neck | none | tilt1deg_gain1pct | 1.54e+03 | 0.209 [0.187-0.223] | 0.396 [0.34-0.46] | 0.282 [0.259-0.312] | -0.07 |
| 12-month template | uniform | neck | none | tilt3deg_gain3pct | 1.48e+04 | 0.0222 [0.0201-0.0232] | 0.0456 [0.0361-0.0507] | 0.0318 [0.0263-0.0351] | -0.03 |
| 12-month template | gradient | neck | none | tilt3deg_gain3pct | 1.53e+03 | 0.206 [0.183-0.221] | 0.383 [0.334-0.453] | 0.279 [0.258-0.308] | -0.07 |
| 12-month template | uniform | neck | homogeneous | tilt0deg_gain0pct | 5.93e-12 | - | - | - | +0.00 |
| 12-month template | gradient | neck | homogeneous | tilt0deg_gain0pct | 1.14e+03 | 0.207 [0.203-0.215] | 0.386 [0.37-0.421] | 0.27 [0.26-0.286] | -0.13 |
| 12-month template | uniform | neck | homogeneous | tilt1deg_gain1pct | 227 | 0.425 [0.367-0.488] | 0.863 [0.771-0.954] | 0.591 [0.557-0.632] | -0.02 |
| 12-month template | gradient | neck | homogeneous | tilt1deg_gain1pct | 1.14e+03 | 0.207 [0.202-0.215] | 0.383 [0.366-0.421] | 0.27 [0.261-0.288] | -0.13 |
| 12-month template | uniform | neck | homogeneous | tilt3deg_gain3pct | 722 | 0.136 [0.125-0.152] | 0.267 [0.248-0.293] | 0.206 [0.19-0.215] | -0.02 |
| 12-month template | gradient | neck | homogeneous | tilt3deg_gain3pct | 1.14e+03 | 0.206 [0.201-0.213] | 0.379 [0.36-0.414] | 0.268 [0.255-0.283] | -0.12 |
| 12-month template | gradient | origin | homogeneous | tilt0deg_gain0pct | 1.14e+03 | 0.207 [0.203-0.215] | 0.386 [0.37-0.421] | 0.27 [0.26-0.286] | -0.13 |
| 12-month template | gradient | origin | homogeneous | tilt1deg_gain1pct | 1.14e+03 | 0.207 [0.203-0.215] | 0.383 [0.368-0.423] | 0.27 [0.261-0.288] | -0.13 |
| 12-month template | gradient | origin | homogeneous | tilt3deg_gain3pct | 1.13e+03 | 0.206 [0.202-0.214] | 0.382 [0.365-0.415] | 0.269 [0.257-0.285] | -0.12 |
| 12-month template | uniform | neck | homogeneous+gradient | tilt0deg_gain0pct | 6.35e-12 | - | - | - | +0.00 |
| 12-month template | gradient | neck | homogeneous+gradient | tilt0deg_gain0pct | 7.06e-13 | - | - | - | +0.00 |
| 12-month template | uniform | neck | homogeneous+gradient | tilt1deg_gain1pct | 224 | 0.39 [0.342-0.436] | 0.808 [0.729-0.878] | 0.508 [0.481-0.547] | -0.02 |
| 12-month template | gradient | neck | homogeneous+gradient | tilt1deg_gain1pct | 23.2 | 3.78 [3.37-4.47] | - | 4.64 [4.06-4.88] | -0.03 |
| 12-month template | uniform | neck | homogeneous+gradient | tilt3deg_gain3pct | 707 | 0.132 [0.12-0.142] | 0.261 [0.239-0.277] | 0.164 [0.151-0.177] | -0.02 |
| 12-month template | gradient | neck | homogeneous+gradient | tilt3deg_gain3pct | 69.6 | 1.29 [1.17-1.39] | 2.56 [2.33-2.71] | 1.49 [1.37-1.59] | -0.03 |
| 12-month template | gradient | origin | homogeneous+gradient | tilt0deg_gain0pct | 5.89e-13 | - | - | - | +0.00 |
| 12-month template | gradient | origin | homogeneous+gradient | tilt1deg_gain1pct | 19.2 | - | - | - | -0.03 |
| 12-month template | gradient | origin | homogeneous+gradient | tilt3deg_gain3pct | 57.2 | 1.6 [1.38-1.73] | 3.03 [2.7-3.21] | 1.91 [1.71-2.03] | -0.03 |

## B'. Exact rigid motion over 60 s (adult)

Slow drift up to 5 deg per axis plus in-band jitter of 0.05 deg RMS per axis (measured: 0.050, 0.050, 0.050 deg), in B0 = 2 nT and G = 5 nT/m, calibration errors [1.0, 0.01] (tilt deg, gain). Peak field change at a sensor (drift included): median 127 pT, maximum 274 pT; this offset moves the sensors' operating point (dynamic range, gain), which is not modelled beyond the calibration errors.

| correction | in-band RMS, exact [fT] | linear prediction [fT] | exact / linear, per channel: median [5th-95th percentile] | largest deviation |
|---|---|---|---|---|
| none | 1387.3 | 1387.7 | 1.0003 [0.9996-1.0013] | 0.5 % |
| homogeneous | 364.7 | 364.4 | 1.0002 [0.9984-1.0024] | 0.9 % |
| homogeneous+gradient | 21.3 | 21.3 | 1.0004 [0.9997-1.0011] | 0.6 % |

## Notes

- A head-mounted array moving rigidly with the head keeps its sensor-to-head geometry, so sustained head displacement changes only the SQUID geometry; the OPM's geometry changes only if the cap slips. The 'known' rows are the ideal limit of movement compensation (continuous head-position tracking for the SQUID, a measured slip for the OPM).
- In a perfectly calibrated array any rigid motion changes the readings by a uniform field plus a symmetric traceless gradient in the head frame, which the 8-term projection removes exactly; in a uniform field the change is uniform and the homogeneous projection removes it exactly; translation in a gradient is uniform, rotation in a gradient is not (tests/test_motion.py, finite rotations included). Calibration errors leave residuals proportional to the field change.
- D compares the OPM after its correction with Neuromag without any projection or SSS (none is charged to it), which is conservative for the OPM. The gradient is defined about x_ref = r0 of G2 (0, 0, 40 mm, head frame), which sets the uniform part of a translated gradient.
- A slipped rigid cap must lift where the head is in its way: sensors that would enter the scalp are moved out along their axes without the 5-mm limit of A-OPM-CLEAR (a flexible cap or sliding holders), so the slip rows include the lifts (count and largest lift per row in the table).
- The static room field is treated as in G2 for both systems (removed exactly by the projection in the projected condition); the same calibration errors would also leave part of it, for both systems, which is not modelled here.
- Not modelled: sensor dynamic range, gain change with the operating field and cross-axis projection beyond the declared calibration errors; head-motion statistics of real children; specific movement-compensation algorithms; field changes from moving magnetic material.
- Field strengths, motion amplitudes, pivot and calibration errors are declared sweeps (A-MOT-*), not measurements of a particular room or device; results scale linearly with rotation x field.
