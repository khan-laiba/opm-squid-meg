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

Dense arrays (G2): farthest-point sampling of the scalp above the brow plane, excluding scalp
points more than 2 mm inside the smooth BEM head surface (ear canals and pinna folds), clearance
as above, then pruning to a 17-mm minimum centre spacing (A-OPM-PACK). The densest feasible
single-axis array has 216 sites (`opm_dense`); the channel-budget control `opm204` takes 204 of
them by farthest-point sampling. Every sensing centre lies outside the BEM head surface (at least
3.0 mm for opm99, 4.9 mm for the dense arrays). A 306-channel single-axis array does not fit.

## 3. Forward models — `opmsquid.forward`, `opmsquid.anatomy`

MNE-Python 1.13.2 BEM forward models (single-compartment inner skull, 0.3 S/m, unless a paper
configuration prescribes otherwise; G2 uses the 3-layer 0.3/0.006/0.3 S/m model). Sources are
taken only among usable vertices, at least 2 mm from the 5120-triangle inner-skull mesh
(A-BEM-DIST): closer vertices carry linear-collocation artefacts (lead-field energy up to ~4000x
that of neighbouring vertices within 0.5 mm, none beyond 2 mm). Cortical sources fixed along the cortical normal (cortical
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

* Anatomy and forward: MNE sample subject, white surface, fixed cortical normals; sources are
  chosen among the 307,858 vertices inside the inner skull and at least 2 mm from its 5120-triangle
  mesh (A-BEM-DIST: closer vertices carry BEM artefacts); 3-layer BEM (5120 triangles per surface,
  0.33/0.0042/0.33 S/m, linear collocation); the sample recording's Vectorview geometry and head
  position; 4-point coil integration (study coils 9014/9024). The matched OPM array (99 sites)
  is the NEW column.
* Descriptors: depth to the nearest node of the 2,562-node BEM scalp; orientation to the normal of
  the nearest inner-skull BEM node, folded to 0-90 deg; the paper's 5 mm x 10 deg bins.
* Sources: 3783 dipoles sampled per bin to the paper's counts (the "dipole traces" cannot be
  rebuilt on new anatomy); 600 nAm PCHIP spike through the digitised Fig. 2(b) waveform. Patches
  grown from each dipole along mesh edges within +/-10 deg until the area first exceeds 20 mm^2
  (2955 of 3783 grow), uniform density with a 622-nAm median total (611-677 nAm; paper 612-678).
* Background: 30,786 random usable vertices (10 %), independent band-limited Gaussian moments,
  each peak-normalised to 10 nAm over a stationary 6 s (generated with 3 s of padding each side),
  one realization shared by all sources and arrays. Two levels (U-HU-bglevel): `as_specified`,
  and `fig6_calibrated`, which multiplies every background moment by 0.43, the ratio of the
  paper's Fig. 6 magnetometer baselines at channels 0631/0711/0741 to ours drawn and digitised
  identically (Appendix F of `docs/literature/hunold2016.md`). The scalar is uncertain: other
  reasonable channel choices give 0.40 (gradiometers 0412/0413/0423) to 0.52 (all-channel
  median), which scales the calibrated bin means by 0.82-1.06. It describes the effective
  background level at the sensors; whether the difference from the text's level comes from the
  normalisation, the filtering, the head position or the anatomy cannot be determined.
* SNR: per sensor type, the channel with the largest noise-free spike; background amplitude
  2 mean|hilbert| over the 1 s before onset (analytic signal of the whole trace, cropped);
  primary numerator the noise-free peak-to-peak spike (Appendix E: 0.87-1.10x the printed Fig. 6
  values depending on digitisation handling, against 0.61-0.73x for the literal noisy peak);
  variants peak, noisy peak and noisy peak-to-peak. Bin means as in the paper; unpaired per-bin
  tests (Student's t if both groups pass Shapiro-Wilk, else rank-sum).
* Comparison with the paper, per bin against the digitised colour classes of Figs 4(a) and 5(a)
  (midpoints). Calibrated, p2p: r = 0.96-0.97 (dipoles), 0.93-0.94 (patches); mean ratio
  0.84-0.97 (0.69-1.03 over the calibration range), bins with paper SNR >= 2.5 at 0.92-0.99x,
  weaker bins 0.80-0.94x; the 2.5 classification agrees in 86-95 % of bins; the GM - MM sign
  agrees in all dipole and patch bins with a paper difference of at least one class. The noisy
  peak-to-peak brackets the paper from above (1.09-1.19x; strong bins 1.00-1.07x, weak bins
  1.19-1.28x), so the paper's low-SNR floor near 1 lies between the two numerators. As
  specified, the maps keep their shape (same r) but reach only 0.36-0.42x the paper.
* Fig. 6 spike check (Table 1 tangential examples; left posterior frontal cortex): the
  magnetometer spike agrees (superficial: ours 6.5 pT median, paper 6.5-7.3 pT) but our
  gradiometer spike is 1.6-1.8x the paper's (214 vs 119-134 pT/m); the gradiometer baseline is
  therefore not an independent confirmation of the calibration, and the head position relative
  to the helmet (not reported) may differ.
* Extension (NEW): white sensor noise (SQUID brochure values; OPM 7/15/30 fT/sqrt(Hz)) added to
  the background; spike, background and noise all filtered 0.5-70 Hz (zero phase); p2p numerator.

## 7. Goldenholz cortical SNR maps (G1C, ADAPT; OPM NEW) — `opmsquid.goldenholz`, `scripts/g1c_goldenholz.py`

Configuration: `configs/goldenholz_reference.toml`.

* Sources: a 10-nAm dipole at each of the 307,858 usable white-surface vertices (cortical normal;
  A-BEM-DIST); patches of geodesic radius 10 and 16 mm (Dijkstra along mesh edges) centred on the
  8,124 usable oct-6 vertices, 50 pAm/mm^2 x vertex area, signed sum (median 435 and 1,149
  vertices, 2.85 and 7.49 cm^2; edge-path distances exceed true geodesics, so the areas are ~91-93 %
  of pi r^2 and patch SNRs ~1 dB below a metrically exact construction; the 16 - 10 difference is
  unaffected).
* Forward: 3-layer BEM, skull 0.006 S/m (probable intended value) and 0.06 S/m (as printed),
  T3 coils, accurate integration; the sample recording's SSP (3 vectors) applied to all gains.
  MEG 2443, marked bad in the recording (baseline RMS 23x the gradiometer median), is excluded:
  N = 102 magnetometers, 203 gradiometers, 305 pooled.
* SNR: Eq. 1 with the 1/N factor, reported per channel set (pooling unstated in the paper).
* Noise, modelled: 1,822 independent cortex-normal sources on a 7-mm Poisson-disk grid (3-D
  distance, usable vertices); s_s^2 from the paper's rule (per-type median of recorded /
  (A A^T)_kk, channel-count weighted). Noise, recorded (ADAPT): per-channel variance of the
  -200-0 ms pre-stimulus baselines of the sample task recording (38,357 samples), SSP, 0.5-100 Hz.
  Lead-field diagnostics: no usable column exceeds ~11x the energy of its neighbours; the largest
  single noise source holds 0.3-0.4 % of the modelled noise power.
* Comparison with the paper (soft; the paper reports no distributions): s_s = 2.64 nAm on our
  grid, 1.78 nAm when normalised to the ~4,000 sources of a 7-mm MNE grid (paper 1.6-1.9 nAm);
  focal SNR median -22.0 dB (pooled, all usable vertices; -22.5 dB at the oct-6 centroids; 5-95 %
  -35.1 to -16.5 dB); 55 % of vertices inside the paper's -29 to -19 dB display range. The deep
  medial regions are darkest (67 % below -29 dB in the cingulate, parahippocampal, entorhinal and
  medial orbitofrontal cortex; cingulate median -31 dB), as in the paper's MEG maps, while the
  medial occipital and paracentral cortex are bright; the insula (lateral but deep) is also low
  (-30 dB). Patch 16 mm minus 10 mm: median 5.4 dB (5.8 dB in the mesial temporal lobe) vs 8.4 dB
  for area scaling of our patches and 10 dB quoted by the paper for the mesial temporal lobe
  (modality unclear, A15); cancellation within the larger folded patches explains the shortfall
  from area scaling. Skull 0.06 vs 0.006 S/m (each with its own calibration): median +0.43 dB
  (pooled), 95th percentile of the absolute change 3.5 dB, largest for near-radial sources.
  Maps use the paper's colour limits.
* Extension (NEW): the matched 99-site OPM array. The brain-noise sources are calibrated on the
  recorded minus empty-room variance (instrument noise is 6 % of the recorded variance for
  magnetometers and 35 % for gradiometers), then intrinsic noise is added explicitly (SQUID
  brochure values; OPM swept 7-30 fT/sqrt(Hz)), 0.5-100 Hz. Brain noise only, the OPM array ties
  both SQUID sensor types (Eq. 1 median +0.1 dB vs magnetometers, -0.04 dB vs gradiometers).
  With intrinsic noise, OPM vs gradiometers is +1.3 to +2.4 dB (the gradiometer noise floor),
  and OPM vs magnetometers +0.2 dB at 7 fT/sqrt(Hz) to -1.0 dB at 30 fT/sqrt(Hz) (crossing near
  12 fT/sqrt(Hz)). Equal channel counts or coverage do not explain these differences (review).

## To be written

G2 noise model (background, environment), endpoints and
uncertainty; G3; G4.
