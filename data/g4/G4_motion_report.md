# G4 extension: head motion and OPM slippage (NEW, bounded secondary analysis)

Static fit was studied first (G2, G3B, G4). This analysis adds what the static maps cannot show, on the G3B arrays and noise conventions (Neuromag at top contact; the refitted dense OPM array; intrinsic + brain noise, the G3B primary condition; 10-nAm cortical-normal dipoles; area-weighted medians over cortical targets). It bounds two mechanisms; it does not establish motion robustness. Configuration: `configs/g4_motion.toml`; code: `scripts/g4_motion.py`, `src/opmsquid/motion.py`.

## A. Sustained displacement: the head in the fixed helmet (Neuromag) vs a slipped cap (dense OPM)

'known': the displaced geometry is known (template and noise of the displaced geometry: ideal movement compensation or a known slip). 'mismatched': the template of the reference geometry applied to the displaced data (noise covariance from those data). A head-mounted array moving with the head has no geometry change; only slips appear for it. A wrong-polarity template counts as -60 dB.

| anatomy | system | displacement | known [dB] | mismatched [dB] | share > 3 dB loss (mismatched) | note |
|---|---|---|---|---|---|---|
| adult | Neuromag combined | down 2 mm | -0.08 | -0.13 | 0.00 | nearest magnetometer 22.0 mm |
| adult | Neuromag combined | x+2 mm | +0.01 | -0.05 | 0.00 | nearest magnetometer 19.5 mm |
| adult | Neuromag combined | x-2 mm | -0.00 | -0.06 | 0.00 | nearest magnetometer 20.1 mm |
| adult | Neuromag combined | y+2 mm | -0.02 | -0.07 | 0.00 | nearest magnetometer 20.6 mm |
| adult | Neuromag combined | y-2 mm | +0.02 | -0.04 | 0.00 | nearest magnetometer 19.7 mm |
| adult | Neuromag combined | down 5 mm | -0.19 | -0.53 | 0.00 | nearest magnetometer 23.0 mm |
| adult | Neuromag combined | x+5 mm | +0.02 | -0.33 | 0.00 | nearest magnetometer 18.2 mm |
| adult | Neuromag combined | x-5 mm | -0.00 | -0.38 | 0.00 | nearest magnetometer 19.0 mm |
| adult | Neuromag combined | y+5 mm | -0.03 | -0.35 | 0.00 | nearest magnetometer 18.9 mm |
| adult | Neuromag combined | y-5 mm | +0.04 | -0.28 | 0.00 | nearest magnetometer 18.7 mm |
| adult | Neuromag combined | down 10 mm | -0.38 | -1.64 | 0.05 | nearest magnetometer 23.2 mm |
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
| adult | Neuromag combined | yaw +10 deg | -0.00 | -1.21 | 0.15 | nearest magnetometer 19.9 mm |
| adult | Neuromag combined | yaw -10 deg | +0.01 | -1.19 | 0.14 | nearest magnetometer 18.2 mm |
| adult | OPM dense | slip x +1 deg | -0.02 | -0.06 | 0.00 | sensors moved 1.9 mm (median); 18 lifted (max 2.5 mm) |
| adult | OPM dense | slip x -1 deg | +0.00 | -0.05 | 0.00 | sensors moved 1.9 mm (median); 30 lifted (max 1.0 mm) |
| adult | OPM dense | slip y +1 deg | -0.01 | -0.06 | 0.00 | sensors moved 1.7 mm (median); 34 lifted (max 2.0 mm) |
| adult | OPM dense | slip y -1 deg | -0.01 | -0.05 | 0.00 | sensors moved 1.7 mm (median); 32 lifted (max 1.5 mm) |
| adult | OPM dense | slip z +1 deg | -0.00 | -0.03 | 0.00 | sensors moved 1.5 mm (median); 13 lifted (max 2.5 mm) |
| adult | OPM dense | slip z -1 deg | -0.00 | -0.02 | 0.00 | sensors moved 1.5 mm (median); 15 lifted (max 1.5 mm) |
| adult | OPM dense | slip x +3 deg | -0.07 | -0.46 | 0.01 | sensors moved 5.6 mm (median); 64 lifted (max 7.5 mm) |
| adult | OPM dense | slip x -3 deg | -0.02 | -0.40 | 0.01 | sensors moved 5.5 mm (median); 74 lifted (max 3.0 mm) |
| adult | OPM dense | slip y +3 deg | -0.06 | -0.42 | 0.01 | sensors moved 5.0 mm (median); 84 lifted (max 5.0 mm) |
| adult | OPM dense | slip y -3 deg | -0.05 | -0.42 | 0.01 | sensors moved 5.0 mm (median); 79 lifted (max 4.0 mm) |
| adult | OPM dense | slip z +3 deg | -0.02 | -0.23 | 0.00 | sensors moved 4.5 mm (median); 49 lifted (max 6.0 mm) |
| adult | OPM dense | slip z -3 deg | -0.02 | -0.19 | 0.00 | sensors moved 4.5 mm (median); 43 lifted (max 4.0 mm) |
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
| 2-year template | OPM dense | slip x +1 deg | -0.02 | -0.05 | 0.00 | sensors moved 1.6 mm (median); 3 lifted (max 1.5 mm) |
| 2-year template | OPM dense | slip x -1 deg | +0.01 | -0.02 | 0.00 | sensors moved 1.6 mm (median); 2 lifted (max 0.5 mm) |
| 2-year template | OPM dense | slip y +1 deg | -0.00 | -0.03 | 0.00 | sensors moved 1.5 mm (median); 6 lifted (max 0.5 mm) |
| 2-year template | OPM dense | slip y -1 deg | -0.00 | -0.03 | 0.00 | sensors moved 1.5 mm (median); 3 lifted (max 0.5 mm) |
| 2-year template | OPM dense | slip z +1 deg | +0.00 | -0.01 | 0.00 | sensors moved 1.3 mm (median); 2 lifted (max 0.5 mm) |
| 2-year template | OPM dense | slip z -1 deg | -0.00 | -0.01 | 0.00 | sensors moved 1.3 mm (median); 2 lifted (max 0.5 mm) |
| 2-year template | OPM dense | slip x +3 deg | -0.06 | -0.38 | 0.02 | sensors moved 4.9 mm (median); 28 lifted (max 4.0 mm) |
| 2-year template | OPM dense | slip x -3 deg | -0.01 | -0.28 | 0.01 | sensors moved 4.7 mm (median); 55 lifted (max 2.0 mm) |
| 2-year template | OPM dense | slip y +3 deg | -0.04 | -0.31 | 0.02 | sensors moved 4.2 mm (median); 49 lifted (max 2.5 mm) |
| 2-year template | OPM dense | slip y -3 deg | -0.04 | -0.32 | 0.02 | sensors moved 4.3 mm (median); 50 lifted (max 2.5 mm) |
| 2-year template | OPM dense | slip z +3 deg | -0.00 | -0.12 | 0.00 | sensors moved 3.9 mm (median); 12 lifted (max 2.0 mm) |
| 2-year template | OPM dense | slip z -3 deg | -0.01 | -0.12 | 0.00 | sensors moved 3.9 mm (median); 14 lifted (max 1.5 mm) |
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
| 12-month template | OPM dense | slip x +1 deg | -0.02 | -0.05 | 0.00 | sensors moved 1.6 mm (median); 2 lifted (max 1.5 mm) |
| 12-month template | OPM dense | slip x -1 deg | +0.01 | -0.02 | 0.00 | sensors moved 1.6 mm (median); 1 lifted (max 0.5 mm) |
| 12-month template | OPM dense | slip y +1 deg | -0.00 | -0.04 | 0.00 | sensors moved 1.5 mm (median); 5 lifted (max 0.5 mm) |
| 12-month template | OPM dense | slip y -1 deg | -0.00 | -0.03 | 0.00 | sensors moved 1.5 mm (median); 3 lifted (max 0.5 mm) |
| 12-month template | OPM dense | slip z +1 deg | -0.00 | -0.01 | 0.00 | sensors moved 1.3 mm (median); 1 lifted (max 0.5 mm) |
| 12-month template | OPM dense | slip z -1 deg | +0.00 | -0.01 | 0.00 | sensors moved 1.3 mm (median); 0 lifted (max 0.0 mm) |
| 12-month template | OPM dense | slip x +3 deg | -0.06 | -0.40 | 0.01 | sensors moved 4.8 mm (median); 32 lifted (max 4.0 mm) |
| 12-month template | OPM dense | slip x -3 deg | -0.01 | -0.31 | 0.02 | sensors moved 4.6 mm (median); 48 lifted (max 2.0 mm) |
| 12-month template | OPM dense | slip y +3 deg | -0.04 | -0.34 | 0.02 | sensors moved 4.2 mm (median); 44 lifted (max 2.5 mm) |
| 12-month template | OPM dense | slip y -3 deg | -0.05 | -0.34 | 0.01 | sensors moved 4.3 mm (median); 46 lifted (max 3.0 mm) |
| 12-month template | OPM dense | slip z +3 deg | -0.01 | -0.13 | 0.00 | sensors moved 3.8 mm (median); 12 lifted (max 1.0 mm) |
| 12-month template | OPM dense | slip z -3 deg | -0.01 | -0.13 | 0.00 | sensors moved 3.8 mm (median); 11 lifted (max 1.0 mm) |

## B. In-band motion of the head-mounted array in a static residual field

The dense OPM array moves rigidly with the head (rotation about a pivot 60 mm below the head origin) in a static field: 1 nT uniform, or a 1 nT/m symmetric traceless gradient (Frobenius norm), random isotropic draws. The artefact covariance is J Sigma J^T for in-band rotation with the given RMS per axis; it is added to the noise after the correction ('none'; 'homogeneous': 3-term homogeneous-field projection; 'homogeneous+gradient': the 8-term projection of G2), applied to signal and noise alike. Calibration errors: RMS tilt of each sensitive axis and RMS gain error, unknown to the analyst. Everything is linear in rotation x field: for a field of B nT the rotation thresholds below divide by B. Neuromag is fixed in the room and has no such term; its static detectability is the comparator for D.

Static reference (no motion): OPM detectability change from the correction alone and D (dense OPM vs Neuromag combined, both intrinsic + brain):

| anatomy | correction | OPM change [dB] | D [dB] |
|---|---|---|---|
| adult | none | +0.00 | +1.00 |
| adult | homogeneous | -0.04 | +0.87 |
| adult | homogeneous+gradient | -0.24 | +0.67 |
| 2-year template | none | +0.00 | +1.84 |
| 2-year template | homogeneous | -0.08 | +1.57 |
| 2-year template | homogeneous+gradient | -0.38 | +1.25 |
| 12-month template | none | +0.00 | +2.04 |
| 12-month template | homogeneous | -0.07 | +1.84 |
| 12-month template | homogeneous+gradient | -0.37 | +1.50 |

In-band artefact per channel (median over channels, then over 32 draws) for 1 deg RMS rotation per axis in the unit field, after each correction [fT]. Then the in-band rotation (deg RMS per axis, in the unit field) at which the median OPM detectability falls by 1 or 3 dB, or D falls to 0 dB, when the artefact is not part of the analyst's noise model (matched filter of the static covariance applied to data that contain it), from the median curve over draws, with the 10th-90th percentiles of the per-draw thresholds (draws that do not reach the level within the tested rotations count as beyond them; '-': not reached up to 5 deg). The draws are common random numbers: the same fields and calibration errors in every correction, pivot and anatomy. Last column: the loss at 5 deg when the artefact is part of the known noise covariance (the optimal filter nulls its at most 3 spatial patterns: the bound for data-driven nulling or regression on measured head motion).

| anatomy | field | pivot | correction | calibration (tilt, gain) | artefact [fT per deg] | 1 dB | 3 dB | D = 0 | oracle loss at 5 deg [dB] |
|---|---|---|---|---|---|---|---|---|---|
| adult | uniform | neck | none | tilt0deg_gain0pct | 1.49e+04 | 0.0168 [0.0129-0.0192] | 0.0314 [0.0255-0.0347] | 0.0158 [0.0123-0.0188] | -0.02 |
| adult | gradient | neck | none | tilt0deg_gain0pct | 1.76e+03 | 0.139 [0.123-0.158] | 0.272 [0.244-0.301] | 0.131 [0.116-0.148] | -0.05 |
| adult | uniform | neck | none | tilt1deg_gain1pct | 1.48e+04 | 0.0169 [0.0129-0.0192] | 0.0316 [0.0255-0.0346] | 0.0158 [0.0123-0.0187] | -0.02 |
| adult | gradient | neck | none | tilt1deg_gain1pct | 1.76e+03 | 0.139 [0.122-0.157] | 0.271 [0.244-0.299] | 0.13 [0.116-0.148] | -0.05 |
| adult | uniform | neck | none | tilt3deg_gain3pct | 1.49e+04 | 0.0164 [0.0127-0.0189] | 0.031 [0.0252-0.0342] | 0.0156 [0.0122-0.0184] | -0.02 |
| adult | gradient | neck | none | tilt3deg_gain3pct | 1.76e+03 | 0.138 [0.12-0.155] | 0.27 [0.24-0.296] | 0.129 [0.115-0.147] | -0.04 |
| adult | uniform | neck | homogeneous | tilt0deg_gain0pct | 7.1e-12 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| adult | gradient | neck | homogeneous | tilt0deg_gain0pct | 1.37e+03 | 0.142 [0.134-0.162] | 0.276 [0.264-0.307] | 0.126 [0.119-0.142] | -0.09 |
| adult | uniform | neck | homogeneous | tilt1deg_gain1pct | 235 | 0.461 [0.419-0.508] | 0.915 [0.853-1.01] | 0.385 [0.355-0.424] | -0.01 |
| adult | gradient | neck | homogeneous | tilt1deg_gain1pct | 1.38e+03 | 0.142 [0.134-0.161] | 0.276 [0.264-0.305] | 0.126 [0.12-0.141] | -0.09 |
| adult | uniform | neck | homogeneous | tilt3deg_gain3pct | 701 | 0.153 [0.136-0.169] | 0.293 [0.267-0.316] | 0.133 [0.122-0.142] | -0.01 |
| adult | gradient | neck | homogeneous | tilt3deg_gain3pct | 1.39e+03 | 0.14 [0.131-0.162] | 0.274 [0.259-0.307] | 0.126 [0.119-0.141] | -0.08 |
| adult | gradient | origin | homogeneous | tilt0deg_gain0pct | 1.37e+03 | 0.142 [0.134-0.162] | 0.276 [0.264-0.307] | 0.126 [0.119-0.142] | -0.09 |
| adult | gradient | origin | homogeneous | tilt1deg_gain1pct | 1.38e+03 | 0.142 [0.134-0.162] | 0.276 [0.264-0.306] | 0.126 [0.12-0.141] | -0.09 |
| adult | gradient | origin | homogeneous | tilt3deg_gain3pct | 1.39e+03 | 0.14 [0.132-0.162] | 0.274 [0.26-0.307] | 0.126 [0.119-0.141] | -0.08 |
| adult | uniform | neck | homogeneous+gradient | tilt0deg_gain0pct | 7.37e-12 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| adult | gradient | neck | homogeneous+gradient | tilt0deg_gain0pct | 9.58e-13 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| adult | uniform | neck | homogeneous+gradient | tilt1deg_gain1pct | 232 | 0.418 [0.375-0.475] | 0.852 [0.784-0.936] | 0.286 [0.269-0.303] | -0.01 |
| adult | gradient | neck | homogeneous+gradient | tilt1deg_gain1pct | 27.8 | 3.17 [2.87-3.57] | - [> 5-> 5] | 2.32 [2.18-2.44] | -0.02 |
| adult | uniform | neck | homogeneous+gradient | tilt3deg_gain3pct | 691 | 0.141 [0.131-0.152] | 0.275 [0.259-0.293] | 0.108 [0.102-0.114] | -0.01 |
| adult | gradient | neck | homogeneous+gradient | tilt3deg_gain3pct | 83.4 | 1.15 [1.06-1.27] | 2.29 [2.12-2.52] | 0.795 [0.729-0.912] | -0.02 |
| adult | gradient | origin | homogeneous+gradient | tilt0deg_gain0pct | 8.59e-13 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| adult | gradient | origin | homogeneous+gradient | tilt1deg_gain1pct | 24 | 3.81 [3.38-4.14] | - [> 5-> 5] | 2.62 [2.48-2.77] | -0.02 |
| adult | gradient | origin | homogeneous+gradient | tilt3deg_gain3pct | 71.2 | 1.3 [1.19-1.42] | 2.57 [2.37-2.77] | 1 [0.896-1.06] | -0.02 |
| 2-year template | uniform | neck | none | tilt0deg_gain0pct | 1.49e+04 | 0.0213 [0.0177-0.0228] | 0.0412 [0.0327-0.0487] | 0.0264 [0.0229-0.0296] | -0.05 |
| 2-year template | gradient | neck | none | tilt0deg_gain0pct | 1.56e+03 | 0.202 [0.173-0.211] | 0.364 [0.322-0.404] | 0.243 [0.227-0.257] | -0.09 |
| 2-year template | uniform | neck | none | tilt1deg_gain1pct | 1.5e+04 | 0.0213 [0.0178-0.0229] | 0.0412 [0.0328-0.0494] | 0.0264 [0.0229-0.0297] | -0.05 |
| 2-year template | gradient | neck | none | tilt1deg_gain1pct | 1.55e+03 | 0.202 [0.173-0.211] | 0.363 [0.321-0.402] | 0.243 [0.227-0.257] | -0.09 |
| 2-year template | uniform | neck | none | tilt3deg_gain3pct | 1.49e+04 | 0.0212 [0.0178-0.0225] | 0.0406 [0.0329-0.0475] | 0.0263 [0.0231-0.0292] | -0.05 |
| 2-year template | gradient | neck | none | tilt3deg_gain3pct | 1.56e+03 | 0.195 [0.173-0.209] | 0.349 [0.322-0.392] | 0.241 [0.225-0.254] | -0.09 |
| 2-year template | uniform | neck | homogeneous | tilt0deg_gain0pct | 6.19e-12 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| 2-year template | gradient | neck | homogeneous | tilt0deg_gain0pct | 1.2e+03 | 0.198 [0.182-0.208] | 0.354 [0.334-0.391] | 0.231 [0.222-0.244] | -0.14 |
| 2-year template | uniform | neck | homogeneous | tilt1deg_gain1pct | 236 | 0.443 [0.393-0.504] | 0.89 [0.812-0.992] | 0.545 [0.521-0.57] | -0.02 |
| 2-year template | gradient | neck | homogeneous | tilt1deg_gain1pct | 1.2e+03 | 0.198 [0.181-0.209] | 0.354 [0.332-0.393] | 0.231 [0.221-0.245] | -0.14 |
| 2-year template | uniform | neck | homogeneous | tilt3deg_gain3pct | 711 | 0.146 [0.132-0.161] | 0.283 [0.261-0.305] | 0.192 [0.166-0.203] | -0.02 |
| 2-year template | gradient | neck | homogeneous | tilt3deg_gain3pct | 1.19e+03 | 0.194 [0.177-0.206] | 0.348 [0.327-0.381] | 0.229 [0.22-0.242] | -0.13 |
| 2-year template | gradient | origin | homogeneous | tilt0deg_gain0pct | 1.2e+03 | 0.198 [0.182-0.208] | 0.354 [0.334-0.391] | 0.231 [0.222-0.244] | -0.14 |
| 2-year template | gradient | origin | homogeneous | tilt1deg_gain1pct | 1.19e+03 | 0.198 [0.181-0.209] | 0.354 [0.332-0.393] | 0.231 [0.221-0.244] | -0.14 |
| 2-year template | gradient | origin | homogeneous | tilt3deg_gain3pct | 1.18e+03 | 0.197 [0.179-0.206] | 0.352 [0.33-0.383] | 0.23 [0.222-0.243] | -0.13 |
| 2-year template | uniform | neck | homogeneous+gradient | tilt0deg_gain0pct | 6.5e-12 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| 2-year template | gradient | neck | homogeneous+gradient | tilt0deg_gain0pct | 8.89e-13 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| 2-year template | uniform | neck | homogeneous+gradient | tilt1deg_gain1pct | 232 | 0.396 [0.367-0.478] | 0.818 [0.771-0.94] | 0.426 [0.382-0.47] | -0.02 |
| 2-year template | gradient | neck | homogeneous+gradient | tilt1deg_gain1pct | 23.3 | 3.9 [3.28-4.56] | - [> 5-> 5] | 3.77 [3.36-4.13] | -0.03 |
| 2-year template | uniform | neck | homogeneous+gradient | tilt3deg_gain3pct | 691 | 0.134 [0.122-0.148] | 0.263 [0.244-0.286] | 0.143 [0.127-0.152] | -0.02 |
| 2-year template | gradient | neck | homogeneous+gradient | tilt3deg_gain3pct | 72.9 | 1.27 [1.2-1.44] | 2.53 [2.39-2.79] | 1.25 [1.18-1.35] | -0.03 |
| 2-year template | gradient | origin | homogeneous+gradient | tilt0deg_gain0pct | 7.47e-13 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| 2-year template | gradient | origin | homogeneous+gradient | tilt1deg_gain1pct | 20.2 | - [4.18-> 5] | - [> 5-> 5] | 4.8 [4.3-> 5] | -0.03 |
| 2-year template | gradient | origin | homogeneous+gradient | tilt3deg_gain3pct | 61.1 | 1.53 [1.39-1.67] | 2.93 [2.73-3.13] | 1.54 [1.39-1.61] | -0.03 |
| 12-month template | uniform | neck | none | tilt0deg_gain0pct | 1.51e+04 | 0.0228 [0.0207-0.0242] | 0.0492 [0.0385-0.053] | 0.0334 [0.0274-0.0372] | -0.04 |
| 12-month template | gradient | neck | none | tilt0deg_gain0pct | 1.54e+03 | 0.21 [0.196-0.223] | 0.4 [0.351-0.463] | 0.286 [0.264-0.311] | -0.09 |
| 12-month template | uniform | neck | none | tilt1deg_gain1pct | 1.51e+04 | 0.0228 [0.0207-0.0242] | 0.0491 [0.0385-0.0529] | 0.0335 [0.0275-0.037] | -0.04 |
| 12-month template | gradient | neck | none | tilt1deg_gain1pct | 1.54e+03 | 0.211 [0.196-0.224] | 0.401 [0.35-0.465] | 0.285 [0.264-0.31] | -0.09 |
| 12-month template | uniform | neck | none | tilt3deg_gain3pct | 1.49e+04 | 0.0224 [0.0205-0.0238] | 0.0467 [0.0376-0.0521] | 0.0327 [0.0272-0.036] | -0.04 |
| 12-month template | gradient | neck | none | tilt3deg_gain3pct | 1.53e+03 | 0.208 [0.18-0.223] | 0.389 [0.33-0.461] | 0.282 [0.259-0.31] | -0.08 |
| 12-month template | uniform | neck | homogeneous | tilt0deg_gain0pct | 4.31e-12 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| 12-month template | gradient | neck | homogeneous | tilt0deg_gain0pct | 1.13e+03 | 0.207 [0.202-0.216] | 0.385 [0.365-0.429] | 0.27 [0.259-0.285] | -0.13 |
| 12-month template | uniform | neck | homogeneous | tilt1deg_gain1pct | 231 | 0.415 [0.373-0.475] | 0.848 [0.78-0.936] | 0.583 [0.545-0.618] | -0.02 |
| 12-month template | gradient | neck | homogeneous | tilt1deg_gain1pct | 1.14e+03 | 0.207 [0.202-0.215] | 0.384 [0.366-0.423] | 0.269 [0.258-0.286] | -0.13 |
| 12-month template | uniform | neck | homogeneous | tilt3deg_gain3pct | 709 | 0.139 [0.123-0.152] | 0.272 [0.246-0.293] | 0.205 [0.182-0.213] | -0.02 |
| 12-month template | gradient | neck | homogeneous | tilt3deg_gain3pct | 1.13e+03 | 0.206 [0.201-0.214] | 0.381 [0.359-0.416] | 0.267 [0.257-0.282] | -0.12 |
| 12-month template | gradient | origin | homogeneous | tilt0deg_gain0pct | 1.13e+03 | 0.207 [0.202-0.216] | 0.385 [0.365-0.429] | 0.27 [0.259-0.285] | -0.13 |
| 12-month template | gradient | origin | homogeneous | tilt1deg_gain1pct | 1.13e+03 | 0.207 [0.203-0.215] | 0.384 [0.366-0.424] | 0.269 [0.258-0.285] | -0.13 |
| 12-month template | gradient | origin | homogeneous | tilt3deg_gain3pct | 1.13e+03 | 0.207 [0.201-0.215] | 0.386 [0.36-0.421] | 0.269 [0.259-0.283] | -0.13 |
| 12-month template | uniform | neck | homogeneous+gradient | tilt0deg_gain0pct | 5.33e-12 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| 12-month template | gradient | neck | homogeneous+gradient | tilt0deg_gain0pct | 7.48e-13 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| 12-month template | uniform | neck | homogeneous+gradient | tilt1deg_gain1pct | 231 | 0.397 [0.336-0.438] | 0.819 [0.718-0.883] | 0.505 [0.432-0.528] | -0.02 |
| 12-month template | gradient | neck | homogeneous+gradient | tilt1deg_gain1pct | 23.1 | 3.89 [3.28-4.8] | - [> 5-> 5] | 4.41 [3.9-> 5] | -0.03 |
| 12-month template | uniform | neck | homogeneous+gradient | tilt3deg_gain3pct | 699 | 0.131 [0.117-0.142] | 0.259 [0.234-0.277] | 0.161 [0.139-0.174] | -0.02 |
| 12-month template | gradient | neck | homogeneous+gradient | tilt3deg_gain3pct | 69.9 | 1.29 [1.15-1.46] | 2.55 [2.29-2.82] | 1.42 [1.3-1.62] | -0.03 |
| 12-month template | gradient | origin | homogeneous+gradient | tilt0deg_gain0pct | 6.9e-13 | - [> 5-> 5] | - [> 5-> 5] | - [> 5-> 5] | +0.00 |
| 12-month template | gradient | origin | homogeneous+gradient | tilt1deg_gain1pct | 19.5 | - [4.13-> 5] | - [> 5-> 5] | - [> 5-> 5] | -0.03 |
| 12-month template | gradient | origin | homogeneous+gradient | tilt3deg_gain3pct | 58.2 | 1.57 [1.42-1.76] | 2.99 [2.77-3.25] | 1.81 [1.63-2.06] | -0.03 |

## B'. Exact rigid motion over 60 s (adult)

Slow drift up to 5 deg per axis plus in-band jitter of 0.05 deg RMS per axis (measured: 0.050, 0.050, 0.050 deg), in B0 = 2 nT and G = 5 nT/m, calibration errors [1.0, 0.01] (tilt deg, gain). Peak field change at a sensor (drift included): median 130 pT, maximum 278 pT; this offset moves the sensors' operating point (dynamic range, gain), which is not modelled beyond the calibration errors.

| correction | in-band RMS, exact [fT] | linear prediction [fT] | exact / linear, per channel: median [5th-95th percentile] | largest deviation |
|---|---|---|---|---|
| none | 1394.9 | 1394.9 | 1.0003 [0.9997-1.0013] | 0.2 % |
| homogeneous | 334.8 | 335.2 | 1.0002 [0.9980-1.0028] | 2.7 % |
| homogeneous+gradient | 24.0 | 24.0 | 1.0003 [0.9998-1.0009] | 0.2 % |

## Notes

- A head-mounted array moving rigidly with the head keeps its sensor-to-head geometry, so sustained head displacement changes only the SQUID geometry; the OPM's geometry changes only if the cap slips. The 'known' rows are the ideal limit of movement compensation (continuous head-position tracking for the SQUID, a measured slip for the OPM).
- In a perfectly calibrated array any rigid motion changes the readings by a uniform field plus a symmetric traceless gradient in the head frame, which the 8-term projection removes exactly; in a uniform field the change is uniform and the homogeneous projection removes it exactly; translation in a gradient is uniform, rotation in a gradient is not (tests/test_motion.py, finite rotations included). Calibration errors leave residuals proportional to the field change.
- D compares the OPM after its correction with Neuromag without any projection or SSS (none is charged to it), which is conservative for the OPM. The gradient is defined about x_ref = r0 of G2 (0, 0, 40 mm, head frame), which sets the uniform part of a translated gradient.
- A slipped rigid cap must lift where the head is in its way: sensors that would enter the scalp are moved out along their axes without the 5-mm limit of A-OPM-CLEAR (a flexible cap or sliding holders), so the slip rows include the lifts (count and largest lift per row in the table).
- The static room field is treated as in G2 for both systems (removed exactly by the projection in the projected condition); the same calibration errors would also leave part of it, for both systems, which is not modelled here.
- Not modelled: sensor dynamic range, gain change with the operating field and cross-axis projection beyond the declared calibration errors; head-motion statistics of real children; specific movement-compensation algorithms; field changes from moving magnetic material.
- Field strengths, motion amplitudes, pivot and calibration errors are declared sweeps (A-MOT-*), not measurements of a particular room or device; results scale linearly with rotation x field.
