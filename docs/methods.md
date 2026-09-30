# Methods (living document)

Status labels: **REPRO** reproduction with the paper's definitions; **ADAPT** methodological
adaptation; **NEW** new study choice. Parameter IDs refer to `docs/provenance_register.md`.

## 1. Analytical spherical benchmark (G1A, REPRO) — `opmsquid.sphere`

A current dipole q at radius r_Q in a spherically symmetric conductor (centre at the origin),
oriented tangentially, is observed by ideal point magnetometers measuring the radial field on a
concentric sphere of radius s. Volume currents do not contribute to B_r, so B_r equals the
radial component of the primary (Biot-Savart) field.

* Peak field: paper Eq. 1 (`bmax_radial`), evaluated in a cancellation-free form.
* Independent numerical check (`numerical_peak_radial`): the radial component of the full
  Sarvas (1987) field maximised over the whole sensor sphere (2-D grid, then Nelder-Mead);
  agreement better than 1e-8 relative (`tests/test_sphere.py`).
* Depth d = h - r_Q below the scalp; OPM on-scalp (s = h), SQUID off-scalp (s = h + 18 mm).
* SNR = peak field / sigma with sigma_OPM = eta sigma_SQUID; equal-SNR depth from Eq. 3 by
  root finding. Only eta enters d_eq; absolute SNR needs sigma_SQUID (R-sigma).
* The sphere centre (r_Q = 0) is silent and excluded; the ratio B_OPM/B_SQUID decreases from
  5.3086 at the brain surface to (113/95)^3 = 1.6829 at the centre, so a crossing exists only
  for 1.6829 < eta < 5.3086.
* Figures (`scripts/g1a_jas_benchmark.py`): Fig. 3 by the validated legacy replica (dotted line at
  the raster-matched 27.53 mm, U-J1); Fig. 4 on the paper's grid r_Q = linspace(0, b, 101)[1:]
  with exact Eq. 3 roots as markers (printed and drawn values in `g1a_benchmark.json`, U-J2);
  absolute SNR with sigma_SQUID = 0.3546 pT recovered from the Fig. 3 axis (R-sigma). d_eq found
  directly from SNR curves is the same for sigma_SQUID x 0.5, 1 and 2, so only eta matters.
* Toy experiment (Fig. 6): 30-nAm tangential target at r_Q = 0.6 b and noise dipoles at 0.4 b and
  0.8 b (caption depths 63/47/31 mm, U-J6), SNR = ratio of peak fields, sensor standoff 0-60 mm.

## 2. Sensor models

### 2.1 Neuromag (SQUID) — `opmsquid.neuromag`

Geometry of the 306 channels (102 locations x {1 magnetometer, 2 planar gradiometers}) from the
MNE sample recording with its measured device-to-head and head-to-MRI transforms. Coil types
set to MRN's T3 variants (3014, 3024; D-T3). MNE's 'accurate' coil integration (8 points per
gradiometer, 16 per magnetometer). Units: T (magnetometers), T/m (gradiometers); never mixed
without noise normalisation.

Measured geometry (sample subject): all coil normals point away from the head; the distance from
each magnetometer pickup coil to the MRI scalp along its inward normal is 25.1-41.4 mm (median
30.9 mm). The TRIUX 18-mm coil-to-helmet-surface spacing (HW-18mm) is therefore not the scalp gap.

### 2.2 OPM — `opmsquid.opm` (NEW)

Single-axis magnetometers with a 10-mm cubic sensing volume (A-OPM-CELL; custom MNE coil 9901,
27 integration points), sensing centre 7 mm from the helmet inner surface (A-OPM-STANDOFF) plus
a declared helmet-to-scalp gap (A-OPM-GAP), sensitive axis along the smoothed scalp normal
(A-OPM-AXIS).

Matched-site array: each Neuromag location is projected along its inward coil normal onto the
MRI scalp; sites below the brow plane (face/neck) are dropped (A-OPM-COVER) and sites whose
sensing volume would collide with the ear pinna or brow are moved outward (A-OPM-CLEAR).
Sample subject: 99 of 102 sites; 3 sites moved out by 3-4 mm; nearest-neighbour spacing
21.7-30.6 mm (median 26.4 mm).

## 3. Forward models — `opmsquid.forward`, `opmsquid.anatomy`

MNE-Python 1.13.2 BEM forward models (single-compartment inner skull, 0.3 S/m, unless a paper
configuration prescribes otherwise). Cortical sources fixed along the cortical normal (cortical
patch statistics), or discrete dipoles at arbitrary points. Content-addressed cache keyed by all
inputs. Checks: discrete and surface forwards agree to < 1e-5; the OPM coil in a realistic BEM
forward behaves as expected (`tests/test_forward.py`).

## 4. Noise — `opmsquid.noise`

One-sided spectra; variance = integral of PSD x |H_composite|^2 over [0, fs/2], where the
composite response is that of the operation actually applied (zero-phase filtering squares
|H|^2). Verified against simulated filtered noise to 2 % (`tests/test_noise.py`).

## 5. Metrics — `opmsquid.metrics`

General metrics (peak-channel SNR, mean-power SNR_dB, known-topography detectability) computed
on noise-normalised quantities; rank-aware whitening from the eigen-decomposition of the
correlation matrix D C D; invariance to channel units verified. Paper-specific metrics are kept
in the modules for each paper.

## 6. Hunold depth-orientation spike maps (G1B, ADAPT; OPM NEW) — `opmsquid.hunold`, `scripts/g1b_hunold.py`

Configuration: `configs/hunold_reference.toml` (every value tagged as printed, chosen or digitised).

* Anatomy and forward: MNE sample subject, all 311,994 white-surface vertices inside the inner
  skull, fixed cortical normals; 3-layer BEM (5120 triangles per surface, 0.33/0.0042/0.33 S/m,
  linear collocation); the sample recording's Vectorview geometry and head position; 4-point
  coil integration (study coils 9014/9024). The matched OPM array (99 sites) is the NEW column.
* Descriptors: depth to the nearest scalp BEM node; orientation to the normal of the nearest
  inner-skull BEM node, folded to 0-90 deg; the paper's 5 mm x 10 deg bins.
* Sources: 3783 dipoles sampled per bin to the paper's counts (the "dipole traces" cannot be
  rebuilt on new anatomy); 600 nAm PCHIP spike through the digitised Fig. 2(b) waveform. Patches
  grown from each dipole along mesh edges within +/-10 deg until the area first exceeds 20 mm^2
  (2895 of 3783 grow; 20.0-25.1 mm^2), uniform density with a 622-nAm median total (99 % of
  totals <= 650 nAm).
* Background: 31,199 random vertices (10 %), independent band-limited Gaussian moments, each
  peak-normalised to 10 nAm over a stationary 6 s (generated with 3 s of padding each side), one
  realization shared by all sources and arrays. Two levels: `as_specified` and `fig6_calibrated`,
  which multiplies every background moment by 0.43, the ratio of the paper's Fig. 6 magnetometer
  baselines to ours drawn and digitised identically (Appendix F of
  `docs/literature/hunold2016.md`; U-HU-bglevel). The gradiometer baselines give 0.39
  independently.
* SNR: per sensor type, the channel with the largest noise-free spike; background amplitude
  2 mean|hilbert| over the 1 s before onset (analytic signal of the whole trace, cropped);
  primary numerator the noise-free peak-to-peak spike (Appendix E), variants peak and noisy peak.
  Bin means as in the paper; unpaired per-bin tests (Student's t if both groups pass Shapiro-Wilk,
  else rank-sum).
* Comparison with the paper: the digitised colour classes of Figs 3-4 (midpoints) per bin:
  Pearson r, mean ratio, and agreement on SNR >= 2.5. Calibrated, p2p: r = 0.96-0.97 (dipoles),
  0.93-0.94 (patches); mean ratio 0.88-0.96; bins with paper SNR >= 2.5 at 0.96-1.04x; 2.5
  classification agrees in 84-93 % of bins; GM - MM sign agrees in 98-100 % of bins with a
  paper difference of at least one class. Weak bins (paper < 2.5) are 0.83-0.92x: the paper's
  maps do not fall below ~1, as expected if its numerator included background at low SNR.
  As specified, the maps keep their shape (same r) but reach only 0.38-0.41x the paper.
* Extension (NEW): white sensor noise (SQUID brochure values; OPM 7/15/30 fT/sqrt(Hz)) added
  to the background, spike, background and noise filtered 0.5-70 Hz (zero phase), p2p numerator.

## 7. Goldenholz cortical SNR maps (G1C, ADAPT; OPM NEW) — `opmsquid.goldenholz`, `scripts/g1c_goldenholz.py`

Configuration: `configs/goldenholz_reference.toml`.

* Sources: a 10-nAm dipole at each of the 311,994 valid white-surface vertices (cortical
  normal); patches of geodesic radius 10 and 16 mm (Dijkstra along mesh edges) centred on the
  8,193 valid oct-6 vertices, 50 pAm/mm^2 x vertex area, signed sum (median 438 and 1,158
  vertices, i.e. ~2.8 and ~7.5 cm^2).
* Forward: 3-layer BEM, skull 0.006 S/m (probable intended value) and 0.06 S/m (as printed),
  T3 coils, accurate integration; the sample recording's SSP (3 vectors) applied to all gains.
  MEG 2443, marked bad in the recording (baseline RMS 23x the gradiometer median), is excluded:
  N = 102 magnetometers, 203 gradiometers, 305 pooled.
* SNR: Eq. 1 with the 1/N factor, reported per channel set (pooling unstated in the paper).
* Noise, modelled: 1,838 independent cortex-normal sources on a 7-mm Poisson-disk grid (3-D
  distance); s_s^2 from the paper's rule (per-type median of recorded / (A A^T)_kk, channel-count
  weighted). Noise, recorded (ADAPT): per-channel variance of the -200-0 ms pre-stimulus
  baselines of the sample task recording (38,357 samples), SSP, 0.5-100 Hz.
* Comparison with the paper (soft; the paper reports no distributions): s_s = 2.63 nAm on our
  grid, 1.78 nAm when normalised to the ~4,000 sources of a 7-mm MNE grid (paper 1.6-1.9 nAm);
  focal SNR median -21.9 dB (pooled; 5-95 % -35.0 to -16.4 dB), 55 % of vertices inside the
  paper's -29 to -19 dB display range, medial lobes lowest (cingulate median -31 dB, insula -30
  dB), matching the paper's dark medial MEG maps. Patch 16 mm minus 10 mm: median 5.4 dB (5.8-5.9
  dB in the mesial temporal lobe) vs 8.2 dB for pure area scaling and 10 dB quoted by the paper
  for the mesial temporal lobe (modality unclear, A15); cancellation within the larger folded
  patches explains the shortfall from area scaling. Skull 0.06 vs 0.006 S/m (each with its own
  calibration): median +0.35 to +0.49 dB, 95th percentile of the absolute change 3.4-4.2 dB.
* Extension (NEW): the matched 99-site OPM array with the same calibrated background plus
  intrinsic noise (OPM 15 fT/sqrt(Hz); SQUID brochure values), 0.5-100 Hz: Eq. 1 OPM vs
  magnetometers median -0.1 dB (OPM better at 47 % of vertices), vs gradiometers +1.7 dB (97 %).
  Eq. 1 averages channel SNRs, so these differences also reflect channel counts and coverage.

## To be written

G2 noise model (background, environment), endpoints and
uncertainty; G3; G4.
