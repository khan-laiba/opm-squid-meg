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

## To be written

G1A figures and toy experiment; G1B Hunold adaptation; G1C Goldenholz adaptation; G2 noise model
(background, environment), endpoints and uncertainty; G3; G4.
