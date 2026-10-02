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
a declared helmet-to-scalp gap (A-OPM-GAP), sensitive axis along the normal of the smooth BEM
head surface averaged within 15 mm (A-OPM-AXIS; within 6.6 deg of the local normal at every site, 95 % of sites within
3.6 deg; the largest deviation is at a lower occipital site of the dense array).
Clearance (A-OPM-CLEAR): every sensing centre is at least 6 mm (standoff - 1 mm) from the MRI
scalp and 4 mm from the BEM head surface, and the 10-mm cell is kept at least 1 mm outside the BEM
head surface at sample points: its 27 integration points and its surface, each face sampled at
11 x 11 points (1-mm spacing, corners and edges included; `opm.cell_volume_points`;
point-to-triangle distance over the triangles around the 10 nearest mesh vertices, the cell
oriented as the forward model builds it in the head frame). Every point of the cell surface lies
within 0.71 mm of a sample and the distance is 1-Lipschitz, so the rule itself guarantees only
0.29 mm for the whole cell (with its centre outside and the head surface one closed surface far
larger than a cell, the cell is then outside it), not 1 mm. The exact cube-to-mesh distance of the
final cells (`opm.exact_cell_clearance`: the closest vertex-face and edge-edge pairs over every
triangle that could be closer, with an intersection test) is at least 1.001 mm (matched) and
1.005 mm (dense) in the adult; over the six G3B anatomies every cell is at least 0.96 mm out, and
7 of their 1,502 cells (scaled controls, 18-month template) lie 0.96-1.00 mm out.
Sites are moved outward along their axis in 0.5-mm steps by at most 5 mm, otherwise they are
infeasible and dropped. v2 checked only the
integration points, which reach +-3.87 mm of the +-5-mm faces: the final review (2026-10-01) found
cell corners up to 0.89 mm inside the head surface at 57 dense and 19 matched sites, so v3 applies
the rule to the whole cell: it moves most sites outward (median 1.0 mm; matched sites by at most
2.5 mm) and drops 7 dense and 2 matched sites, after which the greedy packing places a few dense
sites differently. The rule keeps the forward model's outer boundary out of the cell (see section 3,
A-BEM-SKIN); the 27-point cell average agrees with a 1,000-point one to 1e-3 (`tests/test_opm.py`).

Matched-site array: each Neuromag location is projected along its inward coil normal onto the
MRI scalp; sites below the brow plane (face/neck) or within 20 mm of a preauricular point are
dropped (A-OPM-COVER). Sample subject: 95 of 102 sites; 76 moved out by 0.5-5 mm (v2: 19); nearest-
neighbour spacing 23.1-41.1 mm (median 27.0 mm); sensing centres 6.1-10.3 mm (median 7.0 mm) above
the BEM head surface (Neuromag coils: 25-41 mm from the scalp, median 31).

Dense arrays (G2): farthest-point sampling (15-mm spacing) of the scalp in the same coverage
region, excluding scalp points more than 2 mm inside or 4 mm outside the smooth BEM head surface
(ear canals and pinna folds; the pinna) and within 20 mm of the lower edge of the MRI field of
view, then clearance as above and pruning to a 17-mm minimum centre spacing (A-OPM-PACK). The
field-of-view rule matters because the sample MRI is pitched ~36 deg relative to the head frame:
the flat cap where its head surface is cut lies above the brow plane at the back of the head.
The densest array found has 205 sites (`opm_dense`; greedy, not proven maximal; nearest-neighbour
spacing 17.0-31.2 mm, median 19.0 mm; 164 sites moved out by up to 5 mm, v2: 50; centres 6.1-11.5 mm,
median 7.1 mm, above the BEM head surface); the channel-budget control `opm204` takes 204 of them by
farthest-point sampling. The
dense array reaches about as low at the back as the matched array (lowest scalp point MRI z -68 vs
-67 mm): v2's lower-occipital sites (down to -105 mm, one of them isolated beyond the infeasible
occipito-cervical crease) cannot keep their cells clear within the 5-mm shift. A 306-channel
single-axis array appears infeasible under this rule (best spacing found under the v2 rule
13.7-14.2 mm < 17 mm; the v3 rule only removes sites), though that is not proven. The composition and
placement of every array are pinned by `tests/test_g2.py`.

## 3. Forward models — `opmsquid.forward`, `opmsquid.anatomy`

MNE-Python 1.13.2 BEM forward models (single-compartment inner skull, 0.3 S/m, unless a paper
configuration prescribes otherwise; G2 uses the 3-layer 0.3/0.006/0.3 S/m model). Every 3-layer
model subdivides the 5,120-triangle head surface once (20,480 triangles, same flat geometry;
A-BEM-SKIN, v2): on-scalp sensors sit a few millimetres from it, where the field of the coarse
surface is not converged on this head. On the v1 arrays (cells up to 2.4 mm inside the coarse
surface) refining it changed the field at OPM integration points by a median 0.7-0.8 % of the
array field scale at 3-4 mm (95th percentile 8-12 %), more below 2 mm, and the dense-array
headline by -2.4 % (`results/g2/*_v1_arrays.*`: historic v1-array diagnoses recorded at
271a36d+dirty, kept for the record, not reproducible at this baseline). On the v2 arrays (every integration point
>= 1 mm outside) the same comparison gave a median 0.66-0.74 % at 3-4 mm (95th percentile
7.5-9.8 %), 1.3-1.6 % at 2-2.5 mm and 4-5 % below 1.5 mm, and the dense-array headline changed by
-1.6 % (1.147x coarse, 1.129x refined). On the v3 arrays (the whole cell >= 1 mm outside; no
integration point closer than 2.17 mm) it gives a median 0.61-0.62 % at 3-4 mm (95th percentile
6.0-6.7 %) and 0.7-0.8 % at 2-2.5 mm; the dense-array headline changes by -0.1 % (1.110x coarse,
1.109x refined on the 1,000-target subset; matched 0.996x vs 0.996x), while
Neuromag gains change by 0.08 % (`scripts/study_opm_near_mesh.py`,
`scripts/study_bem_skin_refinement.py`). An exact test on a
3-shell sphere shows MNE's BEM is accurate near a regular surface (95th percentile 0.45 % at 2 mm
with 5,120 triangles, 0.07 % with 20,480; `scripts/study_bem_sphere_accuracy.py`); the sample head
surface has larger triangles (median edge 7.3 mm, 95th percentile 13.5 mm) and a scalp as thin as
3.5-4.7 mm in places. A second refinement is too large for a full BEM solution here. At the v3
array geometry the coarse surface still differs from the refined one by more than 2 % on 18 dense
and 9 matched channels (up to 6 % on a single channel), and in the refined model the
10-mm cell and a point sensor agree within 0.5 % (dense array; 0.5 % matched) in the peak channel
of each of 1,000 targets (median 0.02 %). Assuming, as on the sphere, that the error falls with the square of the mesh size
(an assumption: a second refinement does not fit in memory, so it is not verified on this head),
the refined surface is several times more accurate still. (v1 geometry: 1.207x coarse vs 1.178x refined for the dense array.) Sources are
taken only among usable vertices, inside the inner skull and at least 4 mm from its
5120-triangle mesh (A-BEM-DIST; 91.3 % of the valid vertices). Within 0.5 mm of the mesh, linear
collocation gives lead-field energies up to ~4000x those of neighbouring vertices; refining the
mesh to 20,480 triangles changes Neuromag gains by a median 14 % (90th percentile 67 %) at
2-3 mm, 2.5 % (13 %) at 3-4 mm and 0.8 % (3.6 %) at 4-5 mm, OPM gains by about half as much. The
rule removes the cortex nearest the skull (8.7 %), which makes a superficial OPM advantage
conservative. Cortical sources fixed along the cortical normal (cortical patch statistics), or
discrete dipoles at arbitrary points. Content-addressed cache keyed by all
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
  chosen among the 284,935 usable vertices (inside the inner skull and at least 4 mm from its
  5120-triangle mesh, A-BEM-DIST); 3-layer BEM (5120 triangles per surface as in the paper's
  description, except that v2 refines the head surface to 20,480 triangles in every 3-layer model,
  A-BEM-SKIN, which changes the Neuromag gains by 0.08 %; the JSON configuration field still reads
  5120; 0.33/0.0042/0.33 S/m, linear collocation); the sample recording's Vectorview geometry and head position; 4-point coil
  integration (study coils 9014/9024). The matched OPM array (95 sites) is the NEW column.
* Descriptors: depth to the nearest node of the 2,562-node BEM scalp; orientation to the normal of
  the nearest inner-skull BEM node, folded to 0-90 deg; the paper's 5 mm x 10 deg bins.
* Sources: 3783 dipoles sampled per bin to the paper's counts (the "dipole traces" cannot be
  rebuilt on new anatomy); 600 nAm PCHIP spike through the digitised Fig. 2(b) waveform. Patches
  grown from each dipole along mesh edges within +/-10 deg until the area first exceeds 20 mm^2
  (2928 of 3783 grow), uniform density with a 622-nAm median total (611-722 nAm; paper 612-678).
* Background: 28,494 random usable vertices (10 %), independent band-limited Gaussian moments,
  each peak-normalised to 10 nAm over a stationary 6 s (generated with 3 s of padding each side),
  one realization shared by all sources and arrays. Two levels (U-HU-bglevel): `as_specified`,
  and `fig6_calibrated`, which multiplies every background moment by 0.47: the paper's Fig. 6
  magnetometer baselines at channels 0631/0711/0741 divided by ours drawn and digitised
  identically, averaged over 20 independent background realizations because both the paper's
  baseline and one realization of ours are single draws (one realization alone gives 0.44-0.52,
  5-95 %; Appendix F of `docs/literature/hunold2016.md`). Other channel choices and statistics
  give 0.42 (gradiometers 0412/0413/0423) to 0.54 (all-channel median), which scales the
  calibrated bin means by 0.87-1.11. The scalar describes the effective background level at the
  sensors; whether the difference from the text's level comes from the normalisation, the
  filtering, the head position or the anatomy cannot be determined.
* SNR: per sensor type, the channel with the largest noise-free spike; background amplitude
  2 mean|hilbert| over the 1 s before onset (analytic signal of the whole trace, cropped);
  primary numerator the noise-free peak-to-peak spike (Appendix E: 0.87-1.10x the printed Fig. 6
  values depending on digitisation handling, against 0.61-0.73x for the literal noisy peak);
  variants peak, noisy peak and noisy peak-to-peak. Bin means as in the paper; unpaired per-bin
  tests (Student's t if both groups pass Shapiro-Wilk, else rank-sum).
* Comparison with the paper, per bin against the digitised colour classes of Figs 4(a) and 5(a)
  (midpoints). Calibrated, p2p: r = 0.96-0.97 (dipoles), 0.91 (patches); mean ratio 0.78-0.91
  (0.68-1.01 over the calibration range), bins with paper SNR >= 2.5 at 0.87-0.92x, weaker bins
  0.73-0.90x; the 2.5 classification agrees in 80-93 % of bins; the GM - MM sign agrees in all
  dipole and patch bins with a paper difference of at least one class. The noisy peak-to-peak
  matches the paper on average (0.97-1.04x; strong bins 0.91-0.97x, weak bins 1.01-1.16x;
  classification agreement 86-94 %), so the paper's values lie between the two numerators,
  closer to the noisy one; the primary numerator was fixed from Fig. 6 beforehand and is not
  changed to fit the maps. As specified, the maps keep their shape (same r) but reach only
  0.37-0.43x the paper.
* Fig. 6 spike check (Table 1 tangential examples; left posterior frontal cortex): the
  magnetometer spike agrees (superficial: ours 6.5 pT median, paper 6.5-7.3 pT) but our
  gradiometer spike is 1.6-1.8x the paper's (214 vs 119-134 pT/m); the gradiometer baseline is
  therefore not an independent confirmation of the calibration, and the head position relative
  to the helmet (not reported) may differ.
* Extension (NEW): white sensor noise (SQUID brochure values; OPM 7/15/30 fT/sqrt(Hz)) added to
  the background; spike, background and noise all filtered 0.5-70 Hz (zero phase); p2p numerator.

## 7. Goldenholz cortical SNR maps (G1C, ADAPT; OPM NEW) — `opmsquid.goldenholz`, `scripts/g1c_goldenholz.py`

Configuration: `configs/goldenholz_reference.toml`.

* Sources: a 10-nAm dipole at each of the 284,935 usable white-surface vertices (cortical normal;
  A-BEM-DIST); patches of geodesic radius 10 and 16 mm (Dijkstra along mesh edges) centred on the
  7,661 usable oct-6 vertices, 50 pAm/mm^2 x vertex area, signed sum (median 424 and 1,107
  vertices, 2.81 and 7.29 cm^2; edge-path distances exceed true geodesics, so the areas are ~89-91 %
  of pi r^2 and patch SNRs ~1 dB below a metrically exact construction; the 16 - 10 difference is
  unaffected).
* Forward: 3-layer BEM, skull 0.006 S/m (probable intended value) and 0.06 S/m (as printed),
  T3 coils, accurate integration; the sample recording's SSP (3 vectors) applied to all gains.
  MEG 2443, marked bad in the recording (baseline RMS 23x the gradiometer median), is excluded:
  N = 102 magnetometers, 203 gradiometers, 305 pooled.
* SNR: Eq. 1 with the 1/N factor, reported per channel set (pooling unstated in the paper).
* Noise, modelled: 1,736 independent cortex-normal sources on a 7-mm Poisson-disk grid (3-D
  distance, usable vertices); s_s^2 from the paper's rule (per-type median of recorded /
  (A A^T)_kk, channel-count weighted). Noise, recorded (ADAPT): per-channel variance of the
  -200-0 ms pre-stimulus baselines of the sample task recording (38,357 samples), SSP, 0.5-100 Hz.
  Lead-field diagnostics: no usable column exceeds ~11x the energy of its neighbours; the largest
  single noise source holds 0.3-0.4 % of the modelled noise power.
* Comparison with the paper (soft; the paper reports no distributions): s_s = 2.68 nAm on our
  grid, 1.76 nAm when normalised to the ~4,000 sources of a 7-mm MNE grid (paper 1.6-1.9 nAm);
  focal SNR median -22.1 dB (pooled, all usable vertices; -22.6 dB at the oct-6 centroids; 5-95 %
  -35.3 to -16.7 dB); 56 % of vertices inside the paper's -29 to -19 dB display range. The deep
  medial regions are darkest (67 % below -29 dB in the cingulate, parahippocampal, entorhinal and
  medial orbitofrontal cortex; cingulate median -31 dB), as in the paper's MEG maps, while the
  medial occipital and paracentral cortex are bright; the insula (lateral but deep) is also low
  (-30 dB). Patch 16 mm minus 10 mm: median 5.4 dB (5.7 dB in the mesial temporal lobe) vs 8.3 dB
  for area scaling of our patches (areas 2.81 and 7.29 cm^2) and 10 dB quoted by the paper for the mesial temporal lobe
  (modality unclear, A15); cancellation within the larger folded patches explains the shortfall
  from area scaling. Skull 0.06 vs 0.006 S/m (each with its own calibration): median +0.45 dB
  (pooled), 95th percentile of the absolute change 3.4 dB, largest for near-radial sources.
  Maps use the paper's colour limits.
* Extension (NEW): the matched 95-site OPM array. The brain-noise sources are calibrated on the
  recorded minus empty-room variance (instrument noise is 6 % of the recorded variance for
  magnetometers and 35 % for gradiometers), then intrinsic noise is added explicitly (SQUID
  brochure values; OPM swept 7-30 fT/sqrt(Hz)), 0.5-100 Hz. Brain noise only, the OPM array ties
  both SQUID sensor types (Eq. 1 median +0.1 dB vs magnetometers, -0.02 dB vs gradiometers).
  With intrinsic noise, OPM vs gradiometers is +1.3 to +2.4 dB (the gradiometer noise floor),
  and OPM vs magnetometers +0.2 dB at 7 fT/sqrt(Hz) to -1.0 dB at 30 fT/sqrt(Hz) (crossing near
  12 fT/sqrt(Hz)). Equal channel counts or coverage do not explain these differences (review).

## 8. Realistic adult OPM vs Neuromag comparison (G2, NEW) — `opmsquid.g2`, `scripts/g2_adult_comparison.py`

Configuration: `configs/g2_adult.toml`; assumptions A-G2-*, D-G2-COND in the register. Report:
`results/g2/G2_report.md` (generated from the JSON by `scripts/g2_report.py`).

Design
* Arrays: Neuromag T3 at the sample recording's measured head position (channel sets mag, grad and
  combined, with the full cross-type covariance); OPM arrays `opm_matched` (95 matched sites,
  coverage control), `opm204` (channel-budget control) and `opm_dense` (205 sites, full system),
  all single-axis, 10-mm cell, 7-mm standoff, no extra scalp gap in the primary run (section 2.2).
* Targets: 10-nAm dipoles along the cortical normal at the 7,661 usable oct-6 vertices; geodesic
  patches of 5, 10 and 20 mm radius around each (signed sums; fixed total 10 nAm or fixed density
  0.25 nAm/mm^2). All arrays see the same targets.
* Noise, identical sources for every array, one analysis band (1-40 Hz, zero-phase Butterworth
  order 4; variance = ASD^2 x ENBW of the composite response, 35.1 Hz): intrinsic white noise
  (SQUID brochure values; OPM 15 fT/sqrt(Hz) primary, 7-30 swept); cortical background (1,755
  area-weighted sources on a 7-mm grid, scale fitted once so the median good-gradiometer variance
  equals the measured task-baseline minus empty-room variance in the band; correlated extension
  with lambda = 5 and 10 mm, refitted); room field (8-term external expansion fitted to the
  empty-room recording, through each array's own coil integration points).
* Conditions: intrinsic; intrinsic+brain; intrinsic+brain+env; projected (the 8-dim external
  subspace removed by the same noise-weighted projection for every array, rank n - 8).
* Metrics: peak-channel SNR, mean-power SNR (dB) and known-topography detectability
  sqrt(s^T C^+ s) (rank-aware whitening of the correlation-normalised covariance) with the oracle
  covariance and with plug-in Ledoit-Wolf estimates from 10 and 60 s of data (matched filter from
  the estimate, evaluated under the true covariance). Comparisons: median over targets of
  log2(d_OPM / d_SQUID) with 95 % CIs from a bootstrap over 70 groups, the 68 Desikan-Killiany
  parcels and the two medial-wall ('unknown') labels (neighbouring targets are correlated;
  resampling targets instead gives CIs about 5-11x narrower),
  the share of targets and of parcels where OPM is higher, by depth, distance from the inner
  skull, orientation, lobe and on the cortex. Targets on the medial wall (FreeSurfer 'unknown')
  are kept in every median and reported separately.
* Sensitivity: OPM noise; correlated background; background calibrated on magnetometers; SQUID
  head position (measured; +/-5 mm along each device axis; +/-5 deg pitch; well fitted, 20 mm
  from the nearest magnetometer); OPM scalp gap (0, 3, 6 mm); a joint OPM noise x scalp gap grid;
  plug-in covariance; frequency bands 1-10, 8-30 and 30-80 Hz with a 100-Hz OPM response
  (`scripts/g2_band_sensitivity.py`). Except the joint grid, one factor at a time: the analyses
  show the dependence, they do not bound it. The scalp-gap variants move the primary OPM sites
  outward along their axes (v2; v1 rebuilt the arrays, which added sensors at larger gaps).
* Convergence: background on the 7-mm grid vs every usable vertex; 3-layer vs 1-layer BEM and the
  1-layer BEM refined to 20,480 triangles (1,000 targets); Neuromag 4-point vs accurate coil
  integration; OPM point vs 10-mm cell; oct-6 vs random full-resolution targets; whitening
  tolerance.
* Bridge to G1A: the realistic OPM/magnetometer peak-field ratio vs depth and the equal-SNR depth
  d_eq(eta) (intrinsic noise only, peak-channel SNR) against the sphere benchmark. d_eq comes from
  the depth-binned medians of the ratio, not from per-source crossings.

Results (`results/g2/g2_summary.json`, v3 arrays, run at 19a8fd2: refined head surface and cell
clearance; medians over the 7,661 targets with parcel-bootstrap 95 % CIs; OPM 15 fT/sqrt(Hz);
3-layer BEM)
* Noise validation: the empty-room model (brochure intrinsic noise + fitted room field) gives
  115 fT (magnetometers) and 21.4 fT/cm (gradiometers) against 115 fT and 20.2 fT/cm measured; the
  room field explains 94 % of the magnetometer and 1.6 % of the gradiometer empty-room variance.
  Brain noise (task baseline minus empty room, 1-40 Hz): gradiometers 37.1 fT/cm (calibrated),
  magnetometers 262 fT measured vs 192 fT predicted (0.73x; the independent cortical background
  under-predicts magnetometer noise, likely distant and non-cortical sources). Median brain-noise
  RMS is ~490-510 fT at the OPM sites and 192 fT at the SQUID magnetometers. The model matches one
  scalar of the measured brain noise: per channel the model/measured ratio spans 0.16-1.90
  (gradiometers, 5th-95th percentile, median 0.88) and 0.13-0.78 (magnetometers, median 0.55),
  part of which is measurement noise (one good gradiometer has a negative measured brain variance).
* Not modelled (limitations of every G2-G4 comparison): non-cortical physiological fields
  (cardiac, ocular, muscle), which on-scalp OPMs, being magnetometers, would see at about the
  amplitude of the SQUID magnetometers and which the 8-term projection removes only in their
  uniform-plus-gradient part; the 1/f rise of real OPM noise at low frequencies (white OPM noise is
  assumed; the 30-fT/sqrt(Hz) end of the sweep bounds a uniformly worse sensor, not a coloured
  one); multi-axis OPMs (single-axis sensors: 205 channels against Neuromag's 306, so the
  comparison is not channel-count matched).
* Intrinsic noise only (no brain noise): Neuromag combined beats every OPM array (dense 0.74x
  [0.71-0.76], matched 0.51x; OPM higher for <= 3 % of targets); OPMs beat the gradiometers alone
  (2.2-3.2x). The ratio scales as 1/(OPM noise): break-even at 11.1 (dense) and 7.7 (matched)
  fT/sqrt(Hz) against Neuromag combined.
* With brain noise (intrinsic+brain): matched 95-site OPM 1.00x [0.97-1.02] Neuromag combined
  (higher for 49 % of targets and 39 % of parcels: a tie), OPM 204 1.11x [1.08-1.14], dense
  205-site OPM 1.11x [1.08-1.14] (higher for 90 % of targets and 94 % of parcels); vs gradiometers
  alone 1.10-1.28x, vs magnetometers alone 1.02-1.16x. Adding the room field changes little (0.99,
  1.11, 1.11x). After the external projection (rank n - 8; for Neuromag it acts on all 306 channels
  jointly, so "magnetometers" or "gradiometers" after the projection are subsets of the jointly
  projected data, not standalone sensor-type systems): 0.89x [0.80-0.93], 1.05x, 1.05x
  [0.99-1.10] (dense higher for 57 % of targets; the interval includes 1, so after the projection
  the dense array's advantage is not established); relative to the unprojected room-field condition
  the projection costs the matched array 11 % of its detectability, the dense array 6 % and
  Neuromag 0.5 %.
* Depth (dense vs combined, intrinsic+brain): 1.60x at 10-15 mm, 1.36x at 15-20 mm, 1.21x at
  20-25 mm, 1.12x at 25-30 mm, 1.06x at 30-35 mm and 0.98-1.04x below 35 mm; the matched array is
  1.18x at 10-15 mm and falls to 0.90-0.94x below 40 mm. By distance from the inner skull (dense):
  1.40x at 4-5 mm, 1.07x beyond 10 mm. By lobe (dense): frontal 1.21x, parietal 1.13x, temporal
  1.10x, occipital 1.07x, insula 1.04x, cingulate 1.01x. Without the 558 medial-wall targets
  (7.3 %): dense 1.12x, matched 1.00x. (v1's rise of the dense ratio again below 45 mm, to 1.25x,
  was produced by cells too close to the coarse head surface.)
* Metric dependence: best single channel (peak-channel SNR) dense OPM 0.88x Neuromag (its best
  channel is usually a gradiometer), mean-power SNR 1.06x; detectability 1.11x.
* Estimated covariance (plug-in, Ledoit-Wolf): 10 s of data cost Neuromag combined 15 % and the
  dense OPM 9 % of the oracle detectability (60 s: 3 % and 2 %), so the dense/combined ratio is
  1.19x (10 s) and 1.12x (60 s); matched 1.11x and 1.02x.
* Extended sources: the ratio barely depends on patch size (dense 1.12, 1.12, 1.11x for 5, 10,
  20 mm radius; matched 1.00x), although cancellation reduces the net moment to 0.78, 0.53 and
  0.33 of the scalar moment.
* Sensitivity (dense | matched vs combined, intrinsic+brain): OPM noise 7 -> 30 fT/sqrt(Hz):
  1.24 -> 1.00x | 1.05 -> 0.92x; correlated background (lambda 5, 10 mm): 1.12, 1.14x | 1.00,
  1.01x; background calibrated on magnetometers: 1.11x | 0.99x; head position (+/-5 mm, +/-5 deg,
  well fitted): 1.08-1.14x | 0.98-1.01x; OPM scalp gap 3 and 6 mm (same sites moved out): 1.07 and
  1.02x | 0.98 and 0.95x. Joint OPM noise x scalp gap: dense 0.91x (30 fT/sqrt(Hz), 6 mm) to 1.11x
  (15 fT/sqrt(Hz), 0 mm), matched 0.82-1.00x. Frequency bands (brain scale and room field
  recalibrated per band): dense 1.11x (1-10 Hz), 1.11x (8-30 Hz), 1.09x (30-80 Hz; 1.06x with a
  100-Hz first-order OPM response); matched 0.98-1.00x.
* Convergence: background grid vs every usable vertex changes the median log2 ratios by <= 0.007;
  the 1-layer BEM refined from 5,120 to 20,480 triangles by <= 0.0017; Neuromag 4-point vs
  accurate integration 0.6 % in the gains; OPM 10-mm cell vs point 0.02 % median (95th percentile
  0.1 %; at most 0.5 % in the matched and 0.5 % in the dense array) in the peak-channel field of
  each of the 1,000 convergence targets; oct-6 vs random full-resolution targets <= 0.038;
  whitening tolerance none. Head model: on the convergence subset the dense/combined ratio is
  1.109x (primary), 1.110x with the v1 5,120-triangle head surface and 1.112x with a 1-layer BEM
  (matched 0.996, 0.996, 1.006x); gains differ by 11-12 % between 3 and 1 layers but the ratios
  hardly at all, and the 3- and 1-layer models agree to 0.3 %. In v1 (cells up to 2.4 mm inside the
  coarse head surface) the same comparison gave 1.21x vs 1.13x: the difference was numerical, not
  the skull and scalp.
* Bridge to G1A (intrinsic noise, peak-channel SNR): the realistic OPM/Neuromag magnetometer
  peak-field ratio follows the sphere with the real standoffs (8 mm OPM, 29.8 mm median SQUID);
  the equal-SNR depth at eta = 3 is 27.5 mm (matched) and 30.2 mm (dense) vs 26.8 mm (sphere, real
  standoffs) and 27.7 mm (Jas). OPMs are ahead at every depth for eta <= 2.0 (both arrays) and behind at every depth for eta >= 4.75 (matched) or 5.25 (dense). The idealised
  benchmark's large OPM advantage assumes sensor-noise-limited SNR; with modelled brain noise,
  which on-scalp sensors also see more strongly, it shrinks to about 1.1x overall and vanishes
  for deep sources.

## 9. Epilepsy relevance, adult (G4, NEW) — `opmsquid.ied`, `opmsquid.detection`, `opmsquid.localization`, `scripts/g4_*.py`

Configuration: `configs/g4_epilepsy.toml`; assumptions A-G4-* in the register. The pediatric part
is section 11; head motion and OPM slippage (a bounded secondary extension) are section 12.

Detection design (`scripts/g4_epilepsy_adult.py`)
* Recordings: the G2 noise model in the time domain (A-G4-TS), one realization per 30-s segment
  shared by Neuromag, the matched 95-site OPM and the dense 205-site OPM array (OPM 15
  fT/sqrt(Hz)); 1-40 Hz, 150 Hz after decimation. Identical events (source, strength, morphology,
  time) in every array: 72 locations stratified by depth (10-20, 20-30, 30-45, 45-70 mm) and
  orientation (frontal-heavy: 24 frontal, 14 temporal, 14 parietal, 12 cingulate, 6 insular and 2
  occipital locations), focal dipoles at 10-320 nAm with three spike-wave morphologies, and 10-mm patches;
  1,728 events, one every 2 s.
* Detectors, per array and Neuromag channel set: a known-source/onset oracle (per-trial
  false-positive probability 0.001) and a practical scanner that knows neither time nor source
  (three waveform templates x 716 cortical candidates distinct from the true sources). The
  scanner's templates are the three simulated morphologies, so its absolute sensitivity and
  false-event rates are optimistic; the paired comparison between arrays is less affected. Whitener from
  10 min of baseline null data; thresholds for 1 and 0.2 false events per minute from 20 min of
  independent calibration null data, frozen; held-out null (20 min) gives 0.65-1.15 and 0.1-0.35 false
  events per minute. A hit is an emitted event (a local maximum of the scan statistic, events at
  least the 0.25-s refractory period apart, as counted on the null data) above threshold within
  +/-50 ms of the true spike peak (v3, final review; v2 took the statistic's maximum within the
  window, which also counted peaks that the event rule merges into a nearby higher noise event).
* Statistics: events share locations (18 per depth band, each with 6 strengths x 3 morphologies),
  so the location is the unit. Paired comparisons count, per location, the events detected by only
  one system and test the per-location differences with an exact sign-flip test
  (`detection.sign_flip_p`); the event-level McNemar test is kept for reference only. S50
  intervals come from a bootstrap over locations (1,000 resamples). An S50 outside the tested
  strengths (10-320 nAm) is only bounded (at or below the weakest, beyond the strongest), so every
  resample is kept: the paired ratio of a resample becomes an interval (one-sided, or unrestricted
  when neither system reaches 50 %), and the 95 % interval takes the 2.5th percentile of the lower
  and the 97.5th of the upper bounds, an end at 0 or infinity being open (`detection.censored_interval`;
  the first v3 summaries dropped resamples in which neither system reached 50 %, a review finding).
  The Wilson bands in the figure are event-level and descriptive.

Detection results (v3 arrays, run at 19a8fd2; focal, three morphologies pooled; strength for 50 %
detection, S50, with a bootstrap over locations; practical detector at 1 false event per minute; the
72 locations never lie on the medial wall; p-values uncorrected)
* S50, Neuromag combined vs dense OPM vs matched OPM: 48 vs 37 vs 44 nAm at 10-20 mm, 85 vs 80 vs
  87 nAm at 20-30 mm, 136 vs 135 vs 146 nAm at 30-45 mm, 263 vs 268 vs 320 nAm at 45-70 mm
  (95 % intervals 37-63, 30-55 and 33-61 nAm at 10-20 mm, which overlap; the deepest band's
  upper limits lie beyond the tested 320 nAm). The oracle needs about half the strength (25 vs 16
  vs 23 nAm at 10-20 mm): searching over time and sources costs roughly a factor 2. The intervals
  of separate systems overlap even where the paired comparison is clear; the paired S50 ratio (same
  location resamples for both systems) is the comparison: Neuromag / dense 1.29 [1.03-1.51] at
  10-20 mm, 1.07 [0.95-1.17] and 1.01 [0.93-1.08] in the next two bands, and 0.98 at 45-70 mm, whose
  interval is open at both ends (in 5 % of its resamples one system or neither reaches 50 % within
  the tested strengths).
* Paired, by location (locations with more events detected only by OPM / only by Neuromag
  combined; exact sign-flip p): dense OPM 8/2, 5/1, 4/0 and 0/2 across the four depth bands
  (p = 0.037, 0.19, 0.12, 0.50); matched OPM 8/4, 4/6, 1/3 and 1/9 (p = 0.52, 0.81, 0.62, 0.016:
  worse in the deepest band). With the oracle: dense 14/2, 10/0, 8/2, 8/2 (p = 0.001, 0.002,
  0.22, 0.06); matched 10/4, 3/2, 3/4, 5/5 (p >= 0.48).
* Scoring and run-to-run variability: v2 scored injected spikes by the statistic's maximum near the
  true peak and had other noise draws (an array's channel count changes the random stream). It gave
  dense 16/0, 6/2, 6/1 and 6/1 (p < 0.001, 0.23, 0.28, 0.12; S50 ratio 1.51 [1.25-1.66] at
  10-20 mm) and matched 9/3, 5/4, 4/5 and 0/7 (p = 0.04, 0.85, 0.62, 0.016); an earlier v2 run (one
  dense site fewer) gave dense 13/0, 8/1, 8/2 and 5/2 and matched 9/2, 6/2, 4/5 and 1/8. The dense
  10-20 mm advantage holds in every run, smaller in v3; the deeper location-level results do not
  (20-30 mm: p = 0.03, 0.23, 0.19) and are not claimed; v1's dense advantage in every band (14/0,
  10/2, 12/2, 13/0) is superseded.
* Sensitivity at 1 false event per minute: superficial (10-30 mm) 40-nAm spikes 0.21 (Neuromag),
  0.26 (matched OPM), 0.32 (dense OPM); deep (30-70 mm) 160-nAm spikes 0.34, 0.29, 0.34 (frozen
  1-per-minute thresholds).
* Matched operating points (`scripts/study_g4_matched_rate.py`; the stored simulations
  re-evaluated with every detector's threshold set on the held-out null to 1 false event per
  minute, in-sample for those data; the re-summarised frozen-threshold values reproduce the
  committed ones exactly): dense OPM 9/1, 7/1, 5/0 and 0/2 (p = 0.012, 0.06, 0.06, 0.50; S50 ratio
  1.32 [1.04-1.55] at 10-20 mm); matched OPM 8/4, 5/5, 1/2 and 1/8 (p = 0.44, 1, 1, 0.027). The
  conclusions do not change.
* Consistency with G2 (`scripts/study_g4_vs_g2.py`): at the same 18 locations per band, the G2
  detectability ratio dense / combined (intrinsic+brain+env) has medians 1.43, 1.13, 1.05 and 1.00,
  and the oracle S50 ratio is 1.54, 1.17, 1.06 and 1.14: both favour the dense array most at
  10-20 mm and fall with depth to 30-45 mm; in the deepest band the oracle ratio rises again (8/2
  locations, p = 0.06) while G2's is 1.00. S50 pools the detection curves of 18 locations over a
  factor-2 strength grid, so it need not equal the median ratio.
* Detection agrees with G2: the full OPM system detects superficial spikes at lower strength
  (about a quarter lower at 10-20 mm) and the advantage fades with depth: by location it is
  nominally significant at 10-20 mm for the practical detector (p = 0.04 uncorrected, 0.012 at
  matched held-out rates; ratio 1.29 [1.03-1.51]) and at 10-30 mm for the oracle. A matched-site
  OPM array is not ahead by location in any band (10-20 mm: 8/4, p = 0.52; p = 0.04 in v2 and 0.14
  in the earlier run: not robust) and behind at 45-70 mm (1/9, p = 0.016 uncorrected; the same
  direction in every run).

Localization design (`scripts/g4_localization.py`, bounded)
* 24 locations (2 per depth x orientation stratum), focal dipoles and 10-mm patches at 80 and
  320 nAm, one event per location and condition (96 events), each in 4 s of independent null data
  shared by the arrays; truth as in detection (3-layer BEM, exact geometry).
* Inverse model with bounded mismatch: 1-layer BEM (inner skull) and a coregistration error of
  2 mm and 2 deg, drawn 8 times (translation direction, rotation axis) and shared by all arrays;
  location i uses draw i mod 8 in every condition (v2; v1 cycled draws by event, so each condition
  saw only 2 of them), median displacement at the true source 2.9 mm; noise covariance from 5 min
  of independent null data (Ledoit-Wolf).
* MNE/dSPM (MNE conventions: depth 0.8, SNR 3) on a 5-mm Poisson-disk grid of 3,821 usable
  vertices that excludes the true source vertices; ECD with MNE's `fit_dipole` (same BEM and
  transform, at least 5 mm inside the inner skull) at the spike peak. Errors are measured on the
  MRI through the analyst's (perturbed) transform; the ECD error in the sensor frame is kept as a
  decomposition. Each event also passes through the practical detector (1 false event per
  minute, thresholds from 10 min of independent null data).
* No goodness-of-fit cut: MNE computes GOF on whitened data, where noise adds about one unit per
  channel, so at equal SNR it is higher for arrays with fewer channels (median for detected 320-nAm
  focal events: 72 % matched OPM, 52 % dense OPM, 40 % Neuromag). GOF, chi2/dof and the 95 %
  confidence volume are descriptive.
* Paired OPM-minus-Neuromag comparisons on identical events (one per location): median error
  difference with a bootstrap CI and Wilcoxon signed-rank p; exact McNemar p for joint detection +
  localization within 10 mm; 16 comparisons per OPM array, uncorrected.

Localization results (v3 arrays, run at 19a8fd2, `results/g4/g4_localization_summary.json`;
Neuromag combined, matched OPM, dense OPM)
* Detection in this subset: 320-nAm focal 0.96, 0.79, 0.96; 320-nAm patches 0.67, 0.71, 0.75;
  80-nAm patches 0.08-0.12 (too weak: their "localization" is that of noise).
* ECD, detected events: 320-nAm focal 4.7, 4.2, 4.6 mm on the MRI (3.5, 2.4, 2.4 mm in the sensor
  frame; the coregistration error alone displaces the source by a median 2.9 mm); 320-nAm patches
  7.0, 8.2, 6.2 mm; 80-nAm focal 4.4, 4.8, 5.3 mm. chi2/dof 0.8-1.1: the head-model and
  coregistration mismatch is small against the noise at these SNRs.
* dSPM peak, detected events: 320-nAm focal 12.6, 10.2, 11.4 mm; 320-nAm patches 13.1, 12.9,
  10.8 mm. Support recovery of the 320-nAm patches (share of the N strongest grid sources that are
  members of the patch, N being its number of grid sources; chance ~N/3,821; v2 counted the grid
  sources within the patch radius of its centre instead): medians 0.1, 0.0, 0.2.
* Detected and localized within 10 mm (dSPM | ECD): 320-nAm focal 0.25 | 0.75, 0.38 | 0.67,
  0.42 | 0.71; 320-nAm patches 0.21 | 0.54, 0.29 | 0.50, 0.38 | 0.62.
* Paired (24 locations; 16 comparisons per OPM array, uncorrected p): for 320-nAm patches neither
  OPM array localizes better with dSPM in this run (median error difference 0.0 mm matched,
  p = 0.81; 0.0 mm dense, p = 0.50), and joint detection and localization within 10 mm differs
  little (5/3, p = 0.73; 5/1, p = 0.22). v2 gave -4.8 mm (p = 0.14) and -5.3 mm (p = 0.046), an
  earlier v2 run -4.8 mm (p = 0.003) and -3.7 mm (p = 0.01): the dSPM gain is not robust to the
  noise draws and the v3 arrays. The smallest adult p is a 1.4-mm smaller ECD error of the dense
  array for 320-nAm patches (p = 0.006), which does not survive a correction over the 16
  comparisons (threshold 0.0031). Differences are not excluded.
* With 24 locations the study bounds rather than resolves localization differences. In this model
  the dipole fit is limited by the coregistration error (about 3 mm) for every array; a dSPM
  benefit of on-scalp sensors for extended sources, seen in two earlier runs, is not established.

## 10. Pediatric extension (G3) — `scripts/g3a_jas_size_benchmark.py`, `scripts/g3b_pediatric_helmet.py`

### G3A: head-size benchmark (REPRO; NEW fixed-shell contrast)
* Jas et al. Table 1 heads (h, b): newborn (55, 48), 1 year (70, 62), 8 years (85, 73) and adult
  (95, 80) mm; OPM on the scalp; SQUID on a size-following shell s = h + 18 mm for every head, the
  paper's "size-following, constant-standoff benchmark" (U-J5). Equal-SNR depth from Eq. 3 (exact
  root) and from the paper's grid rule; "normalized d_eq" = depth below the brain surface / b, as
  plotted in Fig. 5B (U-J4); the brain volume where OPM is ahead is reported separately.
* Results at eta = 3 (newborn, 1 year, 8 years, adult): d_eq 30.8, 29.0, 28.1 and 27.7 mm;
  normalized d_eq 49.6, 33.9, 22.0 and 15.8 % (grid rule 49.0, 33.0, 22.0 and 15.0 %; printed:
  50 % newborn, 15 % adult); brain volume with OPM ahead 87.2, 71.2, 52.6 and 40.4 %.
* NEW, idealised: every head concentric in one adult shell (s = 113 mm; gaps 58, 43, 28, 18 mm).
  OPM is ahead throughout the newborn and 1-year brains up to eta = 8.7 and 4.2; for the 8-year
  head normalized d_eq is 49.4 % (87 % of the volume). This concentric sphere is not a helmet
  fit: real heads sit off-centre in a fixed helmet, which G3B models with pediatric anatomy.

### G3B: fixed adult helmet vs head-adaptive OPM (NEW) — `opmsquid.pediatric`, `scripts/g3b_pediatric_helmet.py`
* Anatomies (D-G3-ANAT, A-G3-SCALE). The adult is the G2 subject with the G2 targets and background
  grid. Two size-only controls scale every MRI-frame coordinate of the adult about its MRI origin:
  "school-age size" by 85/95 (Jas et al. Table 1, 8-year vs adult head radius) and "2-year size" by
  the 2-year template's occipitofrontal circumference over the adult's. Their vertices are the
  adult's, so every comparison with the adult is vertex-wise; areas scale with the square of the
  factor and the 4-mm usable-source rule (A-BEM-DIST) is applied on the scaled meshes. The 24-, 18-
  and 12-month infant templates of O'Reilly et al. (2021; `mne.datasets.fetch_infant_template('2yr'
  / '18mo' / '12mo')`, averages of many MRIs of each age, LGPL-2.1, cite O'Reilly et al. 2021 and
  Richards et al. 2016) are used in their native dimensions with their own 3-layer BEMs, dense head
  surfaces, oct-6 source spaces (whose full white surfaces are the full-resolution cortices),
  aparc labels and fiducials (head frames from them). The 2-year template is the primary pediatric
  anatomy (the 2-year size control is scaled to its head circumference); the 18- and 12-month
  templates (added 2026-10-01) show how the results move between averages of one database, not
  between children. They are templates, not individual children: results on them are conditional
  simulations, not population estimates, and anatomical variability is not assessed. The 12-month
  template's skull is 0.25 mm thick at its thinnest (as distributed; 23 of 2,562 inner-skull
  vertices within 1 mm of the outer skull); the 1-layer BEM check does not depend on it. No
  school-aged native anatomy was available (PLAN, inputs).
* Head size from the dense scalp (head frame): occipitofrontal circumference (largest convex-hull
  perimeter of scalp sections parallel to the fiducial plane), breadth, length, vertex height,
  cap volume above the fiducial plane and inter-auricular distance.
* Fixed helmet (A-G3-PLACE). The same Neuromag sensors (coils 3014/3024, 'accurate' integration,
  typical intrinsic noise) for every anatomy. Placements are chosen from the scalp and helmet
  geometry only: "centred" keeps the head frame where the adult's is in the sample recording (the
  ear line in the same place; a smaller head leaves a wide gap at the vertex); "top" raises the
  head (device +z) until the nearest magnetometer coil centre is 20 mm from the scalp (the 18-mm
  Dewar spacing plus 2 mm, G2's 'well fitted' rule; the usual pediatric positioning and the
  primary placement); "back" moves it posterior (device -y, the occiput against the helmet, as for
  a supine child) to the same contact; bounded variants translate the centred head by +-5 mm along
  device x and y, pitch it by +-10 deg or roll it by +-5 deg about the head origin, then raise it
  to the same top contact where there is room (an adult head shifted towards the helmet wall
  stays where it is). After the independent review two variants were added: 'x-centred' shifts the
  centred head along device x until the left and right helmet halves have equal median gaps
  (geometry only), then applies the top contact; 'top-18mm' raises it to true contact at the Dewar
  spacing. A placement is feasible if no
  magnetometer coil centre is within 18 mm of the scalp. Regional gaps use MNE's Vectorview
  selections (left/right frontal, temporal, parietal, occipital).
* Counterfactual helmet (A-G3-COUNTERFACTUAL; a mechanistic control, not a pediatric SQUID
  system): every coil centre scaled by the child/adult head-circumference ratio about the head
  origin of the centred placement (coil sizes, orientations, integration and noise unchanged).
  Around a scaled adult this reproduces the adult's measured fit with every gap scaled by the same
  factor; where a coil would come within 18 mm of the scalp the factor is raised in steps of 0.005.
  It is computed about the centred head and about the laterally centred head ('x-centred').
* Head-adaptive OPM arrays: refitted to each head with the G2 rules unchanged (10-mm cell, 7-mm
  standoff, 17-mm packing, clearance and coverage rules; nothing shrunk): the dense array and the
  matched-site array (the Neuromag sites at the primary placement projected onto the scalp; the
  17-mm packing rule now also applies to it, which removes no site on the adult head, where the
  projected sites are at least 22.9 mm apart). The 204-site channel-budget control has no
  pediatric counterpart: fewer than 204 sites fit on the smaller heads.
* Conventions held fixed (geometry-only comparison first): the G2 analysis band, intrinsic noise
  (OPM 15 fT/sqrt(Hz)), room field, 3-layer BEM conductivities (0.3/0.006/0.3 S/m; the 1-layer BEM
  is the conductivity-free check) and the background moment variance per unit cortical area,
  calibrated once on the adult's Neuromag gradiometers (G2's rule). Background sensitivity: the
  variance x0.5 and x2 (a bounded sensitivity, not an age-specific estimate). The templates'
  averaged white surfaces are smoother than an individual cortex (usable cortex 1,062, 975 and 895
  cm^2 at 24, 18 and 12 months vs 1,878 cm^2 for the adult), which lowers their total background
  power and their patch cancellation.
* Metric (A-G3-METRIC): known-topography detectability d of a 10-nAm cortical-normal dipole with
  the oracle covariance, in dB (20 log10 d); D = dB_OPM - dB_SQUID for each SQUID comparator (102
  magnetometers, 204 gradiometers, combined 306); D_child, D_adult and Delta = D_child - D_adult
  reported together. D does not depend on the moment. Peak-channel and mean-power SNR (dB) are
  secondary metrics. Homologous comparison: vertex-wise for the scaled controls; for the templates
  (no vertex correspondence) per Desikan-Killiany parcel and per declared depth (native mm below
  the scalp) and orientation stratum, from area-weighted medians (target weights = usable cortical
  area of each target's nearest-target cell). Strata or parcels with fewer than 10 targets in
  either anatomy are reported as sparse. Intervals: bootstrap over parcels (one anatomy, or each
  anatomy independently for between-anatomy strata); they contain no between-subject variability.
* Extended sources: fixed-total 10-nAm geodesic patches of 5 and 10 mm (native) around 300 targets
  (the same vertices in the adult and the scaled controls; random targets on each template), for
  the primary placement and the dense OPM array.
* Usefulness (A-G3-USEFUL; operational, not a clinical standard): a source is usable by a system
  when d reaches 5 for a reference moment (20, 50, 100, 200 nAm; maps at 100 nAm); shares of the
  usable cortical area where OPM, SQUID, both or neither are usable, and the moment needed for
  d = 5 by depth.
* The medial wall (FreeSurfer 'unknown': 558 adult targets; 575, 631 and 496 on the 24-, 18- and
  12-month templates) is not cortex and is left out of every G3B summary (it made up most of the
  60-90 mm strata).

G3B results (`results/g3b/g3b_summary.json`, `G3B_report.md`; computed at 19a8fd2 with all six
anatomies and the v3 arrays, A-OPM-CLEAR. Against the v2 pass (cb1b8a9) the adult's D fell from
+0.99 to +0.85 dB with the v3 adult array, while the children's D moved by at most 0.06 dB, so every
Delta rose by 0.12-0.20 dB; the main conclusions did not change, and several within-stratum intervals
now exclude 0, below). Dense OPM vs
Neuromag, intrinsic + brain noise, primary placement unless stated; dB of detectability;
area-weighted medians without the medial wall, parcel-bootstrap 95 % intervals. "Template" alone
means the 2-year template; the 18- and 12-month templates are named.)
* Heads: occipitofrontal circumference 586 mm (adult), 525 (school-age size, x0.895), 495 (2-year
  size, x0.844) and 495 mm (template); usable cortex 1,878, 1,470, 1,290 and 1,062 cm^2. Refitted
  OPM arrays: dense 205, 172, 155 and 151 sites (22.7-25.9 per 100 cm^2 of covered scalp, which
  shrinks from 901 to 584 cm^2), matched 94, 89, 90 and 83 (built from the Neuromag sites at the
  primary placement, which for the adult is 5.5 mm above G2's measured position: one site fewer
  passes the site rules there than G2's 95). Top contact raises the heads by 5.5,
  22, 22 and 28 mm; the median magnetometer-to-scalp gap is then 28.4, 34.1, 39.1 and 38.5 mm
  (centred 29.8, 41.3, 46.9, 47.3 mm). Every placement and the counterfactual helmets are feasible;
  the counterfactual factors are 0.895, 0.844 and 0.904 (template: under the adult's pose it sits
  right of the helmet's midline, centred left/right temporal gaps 48.9/33.3 mm, and its right side
  stops the helmet from shrinking to 0.844). Lateral centring moves the heads by -1.5 (adult), -0.5,
  -0.5 and -6.5 mm (template) along device x; about the laterally centred template the factor is
  0.864 (median gap 27.0 mm, adult 29.8 mm). The 18- and 12-month templates: circumference 491 and
  469 mm, usable cortex 975 and 895 cm^2, dense 157 and 144 sites (26.5 and 25.6 per 100 cm^2 of
  593 and 562 cm^2 covered scalp), matched 80 and 82; top contact raises them by 28 and 33.5 mm
  (median gap 38.9 and 40.4 mm; centred 46.0 and 48.9 mm). They too sit right of the midline
  (centred left/right temporal gaps 47.5/35.5 and 48.1/39.4 mm; lateral centring -5.5 and -2.5 mm),
  so their counterfactual factors are 0.893 and 0.851 (head-circumference ratios 0.838 and 0.801)
  and 0.863 and 0.841 about the laterally centred head.
* Link to G2: at the adult's measured position the dense/combined ratio is 1.111x as G2 reports it
  (unweighted, all targets), 1.124x without the medial wall and +1.09 dB area-weighted; at top
  contact the adult's D is +0.85 dB [+0.61, +1.04] (1.10x).
* D_child, D_adult, Delta (vs Neuromag combined): school-age size +1.45 vs +0.85 dB, Delta +0.59
  [+0.50, +0.71] dB; 2-year size +2.09 vs +0.85, Delta +1.25 [+1.11, +1.37]; 2-year template +1.84
  vs +0.85, Delta +0.88 [+0.55, +1.32] (66 parcels, Delta > 0 in 96 % of their area; the frontal
  poles are sparse). Against the gradiometers alone Delta is +0.80, +1.67 and +1.28 dB, against the
  magnetometers alone +0.56, +1.15 and +0.93 dB. With the external projection +1.10, +1.81 and
  +1.00 dB; for the matched-site OPM array +0.53, +1.30 and +0.93 dB (the gain is not a coverage
  effect of the dense array); for 5- and 10-mm patches +0.69/+0.58, +1.38/+1.29 and +0.69/+0.67 dB.
  Secondary metrics: peak-channel SNR (D_adult -1.30 dB, Neuromag's best channel ahead) Delta
  +0.84, +1.83, +1.72 dB; mean-power SNR (D_adult +0.57 dB) Delta +0.14, +0.53, +0.55 dB.
  18- and 12-month templates: D_child +1.89 and +2.04 dB, Delta +0.99 [+0.58, +1.33] and +1.23
  [+0.73, +1.61] dB (66 and 65 parcels, Delta > 0 in 96 and 97 % of their area); against the
  gradiometers +1.42 and +1.31, the magnetometers +0.73 and +1.21 dB; projected +1.16 and +1.46 dB;
  matched-site array +0.91 and +1.14 dB; 5- and 10-mm patches +1.18/+0.89 and +1.18/+1.11 dB;
  peak-channel SNR +1.55 and +1.69, mean-power SNR +0.67 and +0.86 dB. Across the three templates
  Delta is +0.88, +0.99 and +1.23 dB (24, 18, 12 months; overlapping intervals), but reweighted to
  the adult's depth mix the pooled differences are +0.48, +0.68 and +0.44 dB (below): the
  12-month template's top rank reflects its shallow cortex (the raw ordering also follows head
  circumference, 495, 491 and 469 mm, so depth mix and size are not separated here).
* By depth (template vs adult, combined, native depth strata): Delta +0.26 to +0.67 dB down to 30 mm
  (intervals exclude 0), +0.43 [+0.28, +0.70] at 30-40 mm, +0.49 at 40-50 mm, +1.30 at 50-60 mm and
  +2.55 at 60-90 mm (83 template vs 26 adult targets, 79 of them isthmus cingulate). Within strata the
  template's Delta is largest deep and for radial sources, but those strata hold little of its area
  (2.9 % of it is deeper than 50 mm); most of its pooled gain reflects its shallower cortex (42 % of its area at
  10-20 mm vs 21 % in the adult, area-weighted median depth 21.8 vs 26.2 mm): reweighted to the
  adult's depth mix, the target-level difference of the medians falls from +1.00 to +0.48 dB. For the
  scaled controls the homologous (vertex-wise) Delta is high near the surface and falls with depth
  to 40-50 mm: 2-year size +2.45 [+1.94, +2.67] dB at an adult depth of 10-15 mm, +1.68 at 20-25 mm,
  +0.75 at 30-40 mm, +0.54 at 40-50 mm; school-age size +1.00, +0.82, +0.38 and +0.33 dB; it rises
  again in the deepest strata, which hold few targets (60-90 mm: +1.44 and +1.28 dB). By
  orientation: the scaled controls gain about equally for radial and tangential sources (+0.88 /
  +0.59 dB; +1.39 / +1.33 dB), the template far more for radial sources (0-30 deg +2.48 [+1.82,
  +3.06], 30-60 deg +1.15, 60-90 deg +0.67 dB; radial sources are 16 % of its targets by count). These are
  pooled over depth; at matched depth radial sources still gain more down to 40 mm (+1.00, +0.94,
  +1.03 vs tangential +0.50, +0.35, +0.19 dB at 0-15, 15-25, 25-40 mm) but not deeper (+0.80 vs
  +0.78 dB at 40-90 mm). The 18- and 12-month templates repeat this (`template_depth_checks`):
  within strata Delta is +0.31 to +1.15 dB at 10-50 mm (every interval excludes 0; 12 months at
  10-15 mm +1.15 [+0.41, +1.87]) and +1.23 to +2.45 dB deeper; reweighted to the adult's
  depth mix their pooled difference falls from +1.05 to +0.68 and from +1.20 to +0.44 dB (36 and 46 %
  of their area at 10-20 mm; median depth 23.2 and 20.6 mm); at matched depth radial sources gain
  more than tangential ones down to 40 mm (18 months +1.43, +1.61, +1.28 vs +0.66, +0.55, +0.29 dB;
  12 months +1.73, +0.99, +1.04 vs +0.80, +0.47, +0.18 dB) and only a little more deeper (+0.83 and
  +0.81 vs +0.70 and +0.75 dB at 40-90 mm).
  Regionally the template's Delta is asymmetric (left inferior temporal, entorhinal, fusiform and
  pars orbitalis +2.6 to +3.3 dB; right orbitofrontal -0.7 dB): at top contact under the
  adult's measured pose its left frontal, temporal and parietal gaps are about 8 mm wider than the
  right ones (occipital 1.7 mm); the x-5mm variant, which roughly re-centres it, changes the median
  D by -0.1 dB.
* What drives Delta (each child placement against the adult at the same rule): with the child left
  centred (ear line where the adult's was) Delta is +1.29, +1.85 and +2.06 dB; at top contact +0.59,
  +1.25 and +0.88 dB; laterally centred, then top contact, +0.60, +1.24 and +0.81 dB; at true 18-mm
  contact +0.58, +1.28 and +0.89 dB; back contact +1.03, +1.42 and +1.54 dB. In the counterfactual
  helmet scaled with the head it is -0.08 [-0.11, -0.04], -0.20 [-0.26, -0.16] and +0.43 [+0.19,
  +0.82] dB, and about the laterally centred head -0.08, -0.21 and -0.20 [-0.35, +0.17] dB: the
  template's positive residual came from its lateral offset. For the 18- and 12-month templates:
  centred +2.16 and +2.66, top +0.99 and +1.23, laterally centred then top +0.76 and +1.21, 18-mm
  contact +1.01 and +1.12, back +1.80 and +2.00 dB; counterfactual +0.30 [+0.10, +0.73] and +0.24
  [-0.05, +0.37] dB, about the laterally centred head -0.12 [-0.30, +0.11] and +0.01 [-0.26, +0.25]
  dB. The relative gain of the head-adaptive array in children therefore comes from the fixed
  helmet's fit; with a helmet that fits as the adult's does, it vanishes or reverses slightly.
* Mechanism. Both systems' detectability rises in the smaller heads, by different routes:
  vertex-wise from the adult, the dense OPM gains +1.25 and +1.83 dB in the two scaled controls,
  Neuromag combined +0.65 and +0.54 dB (gradiometers +0.45 and +0.10, magnetometers +0.69 and
  +0.64 dB). The on-scalp OPM sees more signal from a cortex closer in absolute terms (median peak
  field of a 10-nAm dipole 200 fT in the adult, 238, 257 and 271 fT in the children; 255 and 293
  fT at 18 and 12 months) while its brain noise stays near 450-530 fT. In the fixed helmet the
  cortex is farther from the SQUIDs (median target-to-magnetometer distance at 15-20 mm depth 48
  mm adult, 50-57 mm children; OPM 26-27 mm in every head) and, with the background fixed per unit area, the smaller cortex carries
  less background power, so their brain noise falls (magnetometers 203 -> 150, 125 and 120 fT, 110
  and 98 fT at 18 and 12 months; gradiometers 41 -> 30, 23, 23, 21 and 18 fT/cm against 21 fT/cm
  intrinsic) more than their signal (magnetometer peak 70 -> 66, 54, 65, 63 and 67 fT). Absolute
  detectability of a 10-nAm dipole (median dB): dense OPM -1.71 (adult), -0.68, -0.23, +0.92
  (template), +1.08 (18 months) and +1.86 (12 months); Neuromag combined -2.66, -2.15, -2.27,
  -1.19, -1.02 and -0.62. In the counterfactual helmet the SQUID gains at
  least as much as the OPM: the counterfactual Delta depends on the comparator (gradiometers -0.22
  and -0.42 dB, magnetometers -0.03 and -0.12 dB for the scaled controls), which points to a
  SQUID-side component (signal against fixed intrinsic noise, largest for the gradiometers, whose
  brain noise is only about twice their intrinsic noise) rather than to the OPM's fixed standoff
  alone; the noise composition in the counterfactual helmet is not reported.
* Channel count: the children's dense arrays have fewer sites (172, 155, 151). The adult's dense
  array subsampled (farthest-point) to those counts has D +0.63, +0.46 and +0.41 dB (full array
  +0.85), so at an
  equal channel count Delta would be +0.77, +1.56 and +1.51 dB (18- and 12-month templates, 157 and
  144 sites: adult subsampled +0.49 and +0.34, Delta +1.40 and +1.61 dB): the smaller site count of
  the head-adaptive array works against it, and the headline Delta includes that loss.
* Placements: the source-blind variants (+-5 mm, pitch +-10 deg, roll +-5 deg, back contact) give
  D_child 1.41-2.48 (school-age size), 1.64-2.98 (2-year size), 1.75-2.51 (template), 1.80-2.48
  (18 months) and 1.98-2.79 dB (12 months), adult 0.80-1.25 dB; top contact (the primary) ranks in
  the middle of the family. Among the placements, true 18-mm contact gives the lowest D for the
  adult, the school-age size control and the three templates (0.72, 1.32, 1.70, 1.75 and 1.85 dB),
  but not for the 2-year size control (1.98 dB; pitch +10 deg gives 1.64). Regions:
  raising the head (centred -> top) lowers D in every lobe, most in the parietal and frontal
  lobes (template frontal +3.71 -> +2.02, parietal +3.02 -> +1.39 dB); back contact favours the
  SQUID at the occiput (template occipital +1.25 dB vs +2.08 at top) and disfavours it frontally
  (+3.89 vs +2.02 dB).
* Sensitivity (difference of the median D, child minus adult, with the same variant applied to both;
  a check on the medians, not the paired Delta estimator; vs combined): OPM noise 7-30 fT/sqrt(Hz)
  +0.46 to +0.68 (school-age size), +1.09 to +1.31 (2-year size), +0.76 to +1.03 dB (template);
  background variance x0.5 / x2 +0.62 / +0.59, +1.28 / +1.18, +1.04 / +0.94 dB (the variants scale
  the adult too; D_child itself moves by at most 0.06 dB, template 1.86 / 1.79 vs 1.84 dB); 1-layer
  BEM +0.53, +1.17, +0.90 dB. 18- and 12-month templates: OPM noise +0.76 to +1.18 and +0.92 to
  +1.27 dB, background x0.5 / x2 +1.04 / +1.05 and +1.23 / +1.13 dB, 1-layer BEM +0.91 and +1.13 dB.
  At 30 fT/sqrt(Hz) the adult's D is -0.15 dB (Neuromag slightly ahead) and the children's +0.31, +0.94, +0.61, +0.61
  and +0.77 dB.
* Usefulness (d >= 5) at 100 nAm, share of usable cortical area, both / OPM only / SQUID only /
  neither: adult 0.66 / 0.02 / 0.00 / 0.32, school-age size 0.69 / 0.03 / 0.00 / 0.28, 2-year size
  0.69 / 0.05 / 0.00 / 0.26, template 0.75 / 0.05 / 0.00 / 0.20, 18 months 0.76 / 0.05 / 0.00 / 0.19,
  12 months 0.79 / 0.05 / 0.00 / 0.16; at 50 nAm OPM only 0.07, 0.10, 0.14, 0.12, 0.12 and 0.13. Over the
  reference moments (20-200 nAm) the SQUID alone is usable on at most 0.35 % of the area (adult at
  200 nAm; at 100 nAm at most 0.14 %); the cortex neither reaches is mostly deep and
  medial.
* Limits: one adult, three templates of one database and two scaled copies; no between-child
  variability; the templates are averages with smooth cortices; no age-specific background physiology (only the bounded x0.5/
  x2 sensitivity); conductivities held at adult values (the 1-layer BEM bounds this for MEG); the
  intervals contain no between-subject variability. Delta measures a change in relative
  performance under these matching assumptions, not a clinical benefit.

## 11. Epilepsy relevance, pediatric (G4, NEW) — `scripts/g4_epilepsy_pediatric.py`
* The adult detection and localization studies (section 9) run unchanged — configuration
  (`configs/g4_epilepsy.toml`), seeds, strengths, morphologies, detectors, operating points, null
  durations, inverse settings, coregistration error — on a G3B anatomy: Neuromag at the primary
  placement (top contact), the refitted dense and matched OPM arrays, the adult's background
  variance per unit area, room field and intrinsic noise. `g4_epilepsy_adult.Context` holds what
  the studies need (anatomy, arrays, lead fields, background scale); on the adult it reproduces the
  adult results exactly (checked by rerunning the adult through it). Thresholds are calibrated on
  each anatomy's own null data. Locations are stratified by the adult depth and orientation bands;
  a band the anatomy cannot fill gets fewer locations, reported as such.
* Anatomies: the 2-year template (primary), the 18- and 12-month templates, and the two size-only
  controls (scaled adults). Templates are averages and not a population; the comparison with the
  adult mixes head size, anatomy and the fixed-helmet fit, which G3B separates for detectability.

Pediatric G4 results (`results/g4/G4_pediatric_report.md`, `g4_pediatric_comparison.json`; every
anatomy simulated, and the comparison computed, at 19a8fd2 with the v3 arrays; 18 locations per depth
band in every anatomy; p-values uncorrected, over 24 paired detection comparisons (3 comparators x 2
detectors x 4 bands) and 16 localization comparisons per OPM array and anatomy). The children's
Neuromag is at the primary G3B placement (top contact), the adult's at its measured position (the G4
adult study as frozen); top contact would raise the adult's head by only 5.5 mm, against the
children's 22-33.5 mm. Order below: school-age size, 2-year size, 2-year template, 18 months, 12
months.
* Held-out null: 0.4-1.55 false events per minute at the 1-per-minute thresholds (adult 0.65-1.15):
  thresholds calibrated on 20 min and checked on 20 min differ by up to this much, so operating
  points are approximate, and they are lopsided between arrays in either direction (dense OPM vs
  Neuromag combined: 0.4 vs 1.0, 1.25 vs 0.45, 1.3 vs 0.7, 1.1 vs 0.9, 0.85 vs 1.25 per minute).
  At a matched held-out rate of 1 per minute the superficial advantage remains (sensitivity for
  40-nAm spikes at 10-30 mm, dense OPM vs Neuromag combined: 0.31 vs 0.20, 0.31 vs 0.13, 0.43 vs
  0.26, 0.50 vs 0.37, 0.44 vs 0.28; adult 0.33 vs 0.21). With every detector's threshold set on the
  held-out null to 1 false event per minute (`scripts/study_g4_matched_rate.py`,
  `results/g4/G4_matched_rate_report.md`; in-sample for the held-out data) the paired 10-20 mm
  result is unchanged in every anatomy (12-16 locations favour the dense OPM; strength ratio 1.42,
  1.59, 1.40, 1.33 and 1.51; p <= 0.003). Of the deeper differences below, the 2-year size
  control's 45-70 mm (8/1, p = 0.035) and the 12-month template's 30-45 mm (9/1, p = 0.018) remain;
  the 2-year size control's 20-30 and 30-45 mm (9/2, p = 0.06; 7/2, p = 0.11) and the 18-month
  template's 45-70 mm (7/1, p = 0.15) do not. For the matched-site array, 10-20 mm: 6/4, 12/0,
  9/7, 11/1, 11/4 (p = 0.39, 0.0005, 0.66, 0.004, 0.03).
* Strength for 50 % detection, practical detector, Neuromag combined vs dense OPM, 10-20 mm: adult
  48 vs 37, then 50 vs 36, 59 vs 37, 43 vs 30, 36 vs 27 and 43 vs 29 nAm; 45-70 mm adult 263 vs 268,
  then 271 vs 284, 273 vs 243, 240 vs 229, not reached vs 282 and 254 vs 248 nAm. The templates'
  point estimates are lower than the adult's at 10-30 mm, but the intervals overlap (e.g. Neuromag
  at 20-30 mm 71 [61-95], 62 [53-76] and 71 [57-98] vs 85 [68-110] nAm), and the sampled locations'
  median depths differ between anatomies by a few mm (10-20 mm: 13.6-17.8 mm, adult 16.9; 45-70 mm:
  48.0-54.3 mm, adult 51.7): no difference between the anatomies is claimed.
* Paired, dense OPM vs Neuromag combined (locations favouring OPM / Neuromag; paired strength ratio
  Neuromag / OPM): at 10-20 mm every anatomy favours the OPM: adult 8/2 (1.29 [1.03-1.51]), then
  11/1 (1.40 [1.16-1.60]), 16/0 (1.61 [1.33-1.85]), 14/2 (1.44 [1.20-1.65]), 16/0 (1.33
  [1.16-1.62]) and 14/1 (1.46 [1.21-1.84]); p <= 0.04 in each. Deeper, a location-level
  difference appears only for the 2-year size control (20-30 mm 9/2, p = 0.045, 1.19 [0.98-1.46];
  30-45 mm 9/2, p = 0.045, 1.13 [0.99-1.35]; 45-70 mm 9/1, p = 0.02, 1.12 [1.00-open]), the 12-month
  template at 30-45 mm (8/1, p = 0.031, 1.10 [1.03-1.23]) and the 18-month template at 45-70 mm
  (7/1, p = 0.047; Neuromag's S50 is not reached: ratio > 1.13, interval open at both ends); none survives a correction
  over 24 comparisons (smallest p 0.020, threshold 0.0021), and elsewhere p >= 0.12 (ratios
  0.96-1.12). With the oracle every anatomy favours the OPM at 10-20 mm (12/1 to 17/0), five of the
  six at 20-30 mm and two in a deeper band (the 2-year template at 45-70 mm, 11/0, p = 0.001; the
  18-month template at 30-45 and 45-70 mm, 10/3 and 12/2, p = 0.04 and 0.004; the 12-month template
  9/1, 8/2 and 4/5 deeper than 20 mm). The matched-site array favours the OPM at 10-20 mm with the practical detector in the 2-year
  size control (12/0, p = 0.0005) and the 18-month template (12/0, p = 0.0005), not established in
  the 12-month template (11/4, p = 0.051), the school-age control (6/5) or the 2-year template
  (10/5, p = 0.15); like the adult (1/9), the school-age control has a matched-site deficit in its
  deepest band (0/7, p = 0.016).
* So the detectability gains of G3B (Delta +0.59 to +1.25 dB) are not resolved by the practical
  detector with 18 locations per band beyond the superficial band, whose advantage is present in the
  adult and every smaller head and weakest in the adult (strength ratio 1.29 against 1.33-1.61); they
  appear in the oracle's 20-30 mm band and partly deeper.
* Localization (24 locations; Neuromag, matched, dense): ECD errors of detected events are similar
  in every anatomy (320-nAm focal: 4.1-7.9 mm). dSPM, all events, 320-nAm focal: the dense array is
  paired-closer than Neuromag in the 2-year and 18-month templates (-4.7 and -5.3 mm; p = 0.035 and
  0.0033, uncorrected; neither survives the correction below) and the 2-year size control (-1.5 mm, p = 0.01), not in the 12-month template (0.0 mm,
  p = 0.38), the school-age control (0.0 mm, p = 0.25) or the adult (-0.3 mm); 320-nAm patches:
  2-year template -1.2 mm (p = 0.007), 12 months -2.0 mm (p = 0.0061), adult 0.0 mm (p = 0.50).
  ECD, 80-nAm sources, 18 months: -3.0 mm (focal, p = 0.003) and -4.4 mm (patches, p = 0.01).
  Detected and localized within 10 mm (dSPM), 320-nAm focal spikes, Neuromag vs dense: adult 0.25 vs
  0.42; children 0.29 vs 0.33, 0.12 vs 0.29, 0.29 vs 0.54, 0.29 vs 0.58, 0.42 vs 0.58. Within an
  anatomy and OPM array (16 localization comparisons: dSPM and ECD errors and joint
  detection-and-localization success, for focal and patch sources at 80 and 320 nAm; the comparison
  file holds the 8 error comparisons, the localization summaries all 16) four survive a Bonferroni
  correction (p < 0.0031): matched dSPM of 80-nAm patches at 12 months (-33 mm, p = 0.0002; these
  weak patches are mostly not detected, so this compares noise-dominated estimates), matched dSPM of
  320-nAm patches in the school-age control (-10 mm, p = 0.0005), matched dSPM of 320-nAm focal
  events at 18 months (-4.8 mm, p = 0.0019) and dense ECD of 80-nAm focal events at 18 months
  (-3.0 mm, p = 0.0028); across the five children and both arrays (160 comparisons, p < 0.00031)
  only the first. None survives in the 2-year template, the 2-year size control or the adult, and
  which comparisons survive varies between runs (v2 had seven survivors, mostly others). The
  direction is the more robust observation: over all six anatomies 30 of the 192 localization
  comparisons have p < 0.05 (uncorrected; about 10 would be expected by chance if they were
  independent, which they are not), 29 of them in favour of an OPM array; half of the 30 concern
  weak 80-nAm sources, mostly undetected, and 16 are dSPM errors over all events.
* The pediatric epilepsy examples use the same framework as the adult; detection and
  reconstruction claims rest on separate results. Simulated IED-source recovery does not
  identify an epileptogenic zone or establish surgical benefit.

## 12. Head motion and OPM slippage (G4, bounded secondary extension; NEW) — `opmsquid.motion`, `scripts/g4_motion.py`
Static fit is studied first (sections 8-11). GOAL G4 asks for a separate, time-varying analysis of
motion: a head-mounted array keeps its sensor-to-head geometry while its relationship to the residual
room field changes, and a cap can slip. This section bounds both mechanisms on the G3B arrays and
noise conventions (Neuromag at top contact; the refitted dense OPM array; intrinsic + brain noise;
10-nAm cortical-normal dipoles; area-weighted medians over cortical targets) for the adult and the
24- and 12-month templates. It does not establish motion robustness.
* A. Sustained displacement (A-MOT-GEOM). The head moves inside the fixed Neuromag helmet by 2, 5
  and 10 mm (down, +-x, +-y) or rotates by +-5 and +-10 deg (pitch, roll, yaw about the head
  origin); a position that brings a magnetometer coil centre within 18 mm of the scalp is
  infeasible (the head cannot move up from top contact). A head-mounted array moving with the head
  has no geometry change; its cap can slip, modelled as a rigid rotation of the array about the
  head origin by +-1 and +-3 deg about x, y and z, sensors that would enter the scalp being moved
  out along their axes to the A-OPM-CLEAR clearances without its 5-mm limit on the extra shift (a
  flexible cap or sliding holders can lift, not sink; the number of lifted sensors and the largest
  lift are reported per case). For each case,
  'known': the detectability with the displaced geometry (template and noise of the displaced
  geometry; the ideal limit of movement compensation or of a measured slip); 'mismatched': the
  output SNR (h^T s_k) / sqrt(h^T C_k h) of the matched filter h = C_k^+ s_0 built with the
  template s_0 of the reference geometry and the noise covariance C_k of the displaced data
  (`motion.mismatched_detectability`; a wrong-polarity template counts as -60 dB).
* B. In-band motion in a static residual field. A rigid array moving with the head reads
  n_i(t) . B(p_i(t)) at every integration point, with B(x) = B0 + G (x - x_ref) static in the room
  (`motion.readings`, exact for any rotation; `motion.jacobian`, its linearisation in the rotation
  vector about a pivot and the translation). With perfect calibration any rigid motion changes the
  readings by those of a uniform field plus the symmetric traceless gradient R^T G R about x_ref in
  the head frame, so the 8-term projection of G2 removes the change exactly, for any rotation; in a
  uniform field the change is uniform and the homogeneous-field projection removes it exactly;
  translation in G is uniform, rotation in G is not. `tests/test_motion.py` checks these identities
  (finite rotations up to 1 rad included), the Jacobian against finite differences (with a pivot and
  calibration errors), the basis against `environment.external_basis`, the two metrics' bounds and
  the threshold interpolation. Calibration errors (A-MOT-CAL: sensitive-axis tilt 0, 1, 3 deg RMS
  with gain errors 0, 1, 3 % RMS, unknown to the analyst) leave residuals proportional to the field
  change. Head rotation about a pivot 60 mm below the head origin (A-MOT-PIVOT; the head origin as a
  sensitivity), independent per axis with in-band RMS theta; fields of unit strength (A-MOT-FIELD:
  1 nT uniform, or a 1 nT/m symmetric traceless gradient with unit Frobenius norm, about x_ref = r0
  of G2, 0/0/40 mm in the head frame; 32 isotropic draws used as common random numbers: the same
  fields and calibration errors in every correction, pivot and anatomy). The artefact covariance
  theta^2 P J J^T P^T
  after the correction P (none, the 3-term homogeneous or the 8-term projection, applied to signal
  and noise alike) is evaluated two ways (A-MOT-METRIC): outside the analyst's noise model (the
  matched filter of the static covariance applied to data that contain the artefact), and inside
  it (the oracle covariance: a rigid array's linear motion artefact spans at most three spatial
  patterns per field, which the optimal filter nulls; the bound for data-driven nulling or for
  regression on measured head motion). The rotation at which the median OPM detectability falls by
  1 or 3 dB, or to Neuromag's static detectability (D = 0), is interpolated on a logarithmic grid
  (0.0001-5 deg RMS) on the median curve over draws, with the 10th-90th percentiles of the per-draw
  thresholds; everything is linear in rotation x field, so for a field of B nT the thresholds divide
  by B. The SQUIDs are fixed in the room and have no such term; D charges no projection or SSS to
  Neuromag, which is conservative for the OPM.
* B'. A time-domain check with exact rigid motion over 60 s at 300 Hz (adult): slow rotation drift
  (below 0.1 Hz, up to 5 deg per axis) plus in-band jitter (0.05 deg RMS per axis after the G2
  analysis filter), B0 = 2 nT and G = 5 nT/m, calibration errors 1 deg and 1 %; the in-band
  residual after each correction against the linear prediction from the measured in-band motion
  covariance, and the peak field change at a sensor (the shift of its operating point).
* Not modelled: sensor dynamic range, gain changes with the operating field and cross-axis
  projection beyond the declared calibration errors; head-motion statistics of real children;
  specific movement-compensation algorithms; moving magnetic material. The static room field is
  treated as in G2 for both systems; the same calibration errors would leave part of it too.

Motion results (`results/g4/g4_motion_summary.json`, `G4_motion_report.md`, `Figure_G4_motion.png`;
computed at 19a8fd2; medians over cortical targets of the adult and the 24- and 12-month templates;
dB of detectability, dense OPM or Neuromag combined, intrinsic + brain noise)
* A. Neuromag, head displaced in the fixed helmet, template of the reference position: 2 mm costs
  0.04-0.14 dB, 5 mm 0.28-0.56 dB, 10 mm down 1.65-1.72 dB (5-8 % of the cortex losing more than 3
  dB) and 10 mm sideways or forward/back 1.37-1.46 dB (templates only: at top contact the adult's
  head has no room); 5-deg rotations 0.22-0.84 dB, 10-deg rotations 0.91-2.31 dB (pitch and roll of
  the templates 1.90-2.31 dB, with 34-41 % of their cortex losing more than 3 dB). With the displaced
  geometry known (ideal movement compensation) the loss is at most 0.39 dB (10 mm away from the
  helmet top). 18 of the 81 displacements are infeasible at top contact (translations towards the
  helmet wall and some pitch and roll rotations; no upward displacement was tested). Dense OPM cap
  slipped, template of the reference geometry: 1 deg (median sensor shift 1.3-1.9 mm) costs
  0.01-0.09 dB, 3 deg (3.8-5.7 mm) 0.12-0.51 dB; with the slip known at most 0.08 dB. A 3-deg slip
  lifts 11-91 sensors (largest lift 4.5 mm, adult, about x; v2 up to 24 mm), most on the adult's less
  spherical head.
  Per mm of sensor-to-head displacement the uncompensated losses of the two systems are similar; the
  head-mounted array is unaffected by head displacement itself, which costs the fixed helmet up to
  1.7 dB at 10 mm and 2.3 dB at 10 deg without compensation.
* B. In-band rotation, artefact outside the noise model; thresholds in deg RMS per axis for a unit
  field (divide by the residual field in nT or nT/m), 1-dB loss on the median curve over 32 draws,
  ranges over the three anatomies (the 10th-90th percentiles of the per-draw thresholds, a draw
  that does not reach the level within the tested rotations counting as beyond them, lie 2-17 %
  below and 4-23 % above them): no correction 0.021-0.023 deg (uniform) and 0.16-0.21 deg (gradient); homogeneous
  projection: the uniform term is removed exactly with perfect calibration, 0.42-0.44 deg with 1-deg/
  1-% calibration errors and 0.14-0.15 deg with 3 deg/3 %, while the gradient term is not removed
  (0.15-0.21 deg, as without correction); 8-term projection: uniform 0.38-0.40 and 0.13 deg,
  gradient 3.0-3.9 and 1.1-1.3 deg (1 deg/1 % and 3 deg/3 %). The OPM falls to Neuromag's static
  detectability (D = 0) at 0.46-1.5 times these rotations, earliest for the adult, whose static D
  after the projections is smallest (+0.32 to +0.65 dB, against +1.25 to +1.84 dB for the
  templates). With the head origin instead of the neck as pivot the thresholds after the homogeneous
  projection change by less than 2 % (not at all with perfect calibration: the translation the neck
  pivot adds is uniform and removed),
  and after the 8-term projection they rise by 12-22 % (or beyond 5 deg): with calibration errors the
  projection leaks part of the translation term. With the artefact part of the noise
  model (oracle), the loss stays below 0.17 dB up to 5 deg in every case. Artefact per channel for 1
  deg RMS in the unit field: 15 pT without correction (uniform), 0.23-0.24 pT after the homogeneous
  projection with 1 deg/1 % errors, 23-28 fT for the gradient term after the 8-term projection. For
  example, with 1-deg/1-% calibration and the 8-term projection the dense array loses 1 dB at about
  0.4 deg RMS of in-band head rotation in a 1-nT residual field and 0.04 deg in 10 nT; without any
  correction at 0.02 deg in 1 nT.
* B'. Exact rigid motion (adult; 5-deg drift, 0.05-deg in-band jitter, 2 nT and 5 nT/m, 1 deg/1 %):
  per channel the exact in-band artefact is within 0.25 % of the linear prediction for 90 % of the
  channels (largest deviation 0.9 %; medians 1,392, 365 and 22 fT RMS after no, the homogeneous and
  the 8-term correction, against 89 fT of intrinsic OPM noise in the band); the drift moves each
  sensor's operating point by up to 274 pT (median 125 pT).
* Reading. A head-mounted array removes the geometry error that head motion causes in a fixed
  helmet, and converts motion into a field artefact whose size scales with the residual field and
  whose removal depends on calibration (or on modelling it from the data or from measured motion).
  These are bounds on two mechanisms in declared conditions, not an estimate of motion robustness
  in children.
