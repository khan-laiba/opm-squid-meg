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
head surface averaged within 15 mm (A-OPM-AXIS; within 7.1 deg of the local normal at every site, 95 % of sites within
3.5 deg; the largest deviation is at a midline occipital site of the dense array).
Clearance (A-OPM-CLEAR): every sensing centre is at least 6 mm (standoff - 1 mm) from the MRI
scalp and 4 mm from the BEM head surface, and every integration point of the 10-mm cell is at
least 1 mm outside the BEM head surface (exact point-to-triangle distance with the inside test,
the cell oriented as the forward model builds it in the head frame; the first v2 check built it in
the MRI frame, with its in-plane axes rotated by up to ~100 deg, a fault found by the second
review); sites are moved outward along their axis by at most 5 mm, otherwise they are infeasible
and dropped. The cell rule (v2) keeps the
forward model's outer boundary out of the cell (see section 3, A-BEM-SKIN).

Matched-site array: each Neuromag location is projected along its inward coil normal onto the
MRI scalp; sites below the brow plane (face/neck) or within 20 mm of a preauricular point are
dropped (A-OPM-COVER). Sample subject: 97 of 102 sites; 19 moved out by 0.5-5 mm (pinnae, brow,
occiput); nearest-neighbour spacing 22.9-30.3 mm (median 26.6 mm); sensing centres 5.4-10.3 mm
(median 6.3 mm) above the BEM head surface (Neuromag coils: 25-41 mm from the scalp, median 31).

Dense arrays (G2): farthest-point sampling (15-mm spacing) of the scalp in the same coverage
region, excluding scalp points more than 2 mm inside or 4 mm outside the smooth BEM head surface
(ear canals and pinna folds; the pinna) and within 20 mm of the lower edge of the MRI field of
view, then clearance as above and pruning to a 17-mm minimum centre spacing (A-OPM-PACK). The
field-of-view rule matters because the sample MRI is pitched ~36 deg relative to the head frame:
the flat cap where its head surface is cut lies above the brow plane at the back of the head.
The densest array found has 212 sites (`opm_dense`; greedy, not proven maximal; nearest-neighbour
spacing median 18.7 mm; 50 sites moved out by up to 5 mm; centres 5.3-11.5 mm above the BEM head
surface); the channel-budget control `opm204` takes 204 of them by farthest-point sampling. The
dense array reaches lower at the back than the matched array (lowest scalp point MRI z -105 vs
-67 mm), with one isolated lower-occipital site (nearest neighbour 43 mm) beyond the infeasible
occipito-cervical crease. A 306-channel single-axis array appears infeasible under this rule
(best spacing found 13.7-14.2 mm < 17 mm), though that is not proven. The composition and
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
>= 1 mm outside) the same comparison gives a median 0.66-0.74 % at 3-4 mm (95th percentile
7.5-9.8 %), 1.3-1.6 % at 2-2.5 mm and 4-5 % below 1.5 mm; the dense-array headline changes by
-1.6 % (1.147x coarse, 1.129x refined on the 1,000-target subset; matched 1.008x vs 1.004x), while
Neuromag gains change by 0.08 % (`scripts/study_opm_near_mesh.py`,
`scripts/study_bem_skin_refinement.py`). An exact test on a
3-shell sphere shows MNE's BEM is accurate near a regular surface (95th percentile 0.45 % at 2 mm
with 5,120 triangles, 0.07 % with 20,480; `scripts/study_bem_sphere_accuracy.py`); the sample head
surface has larger triangles (median edge 7.3 mm, 95th percentile 13.5 mm) and a scalp as thin as
3.5-4.7 mm in places. A second refinement is too large for a full BEM solution here. At the v2
array geometry the coarse surface still differs from the refined one by more than 2 % on 50 dense
and 26 matched channels (up to 22 % on a single channel), and in the refined model the
10-mm cell and a point sensor agree within 2.8 % (dense array; 1.0 % matched) in the peak channel
of each of 1,000 targets (median 0.03 %). Assuming, as on the sphere, that the error falls with the square of the mesh size
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
  integration (study coils 9014/9024). The matched OPM array (97 sites) is the NEW column.
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
  (-30 dB). Patch 16 mm minus 10 mm: median 5.4 dB (5.7 dB in the mesial temporal lobe) vs 8.5 dB
  for area scaling of our patches and 10 dB quoted by the paper for the mesial temporal lobe
  (modality unclear, A15); cancellation within the larger folded patches explains the shortfall
  from area scaling. Skull 0.06 vs 0.006 S/m (each with its own calibration): median +0.45 dB
  (pooled), 95th percentile of the absolute change 3.4 dB, largest for near-radial sources.
  Maps use the paper's colour limits.
* Extension (NEW): the matched 97-site OPM array. The brain-noise sources are calibrated on the
  recorded minus empty-room variance (instrument noise is 6 % of the recorded variance for
  magnetometers and 35 % for gradiometers), then intrinsic noise is added explicitly (SQUID
  brochure values; OPM swept 7-30 fT/sqrt(Hz)), 0.5-100 Hz. Brain noise only, the OPM array ties
  both SQUID sensor types (Eq. 1 median +0.1 dB vs magnetometers, -0.03 dB vs gradiometers).
  With intrinsic noise, OPM vs gradiometers is +1.3 to +2.4 dB (the gradiometer noise floor),
  and OPM vs magnetometers +0.2 dB at 7 fT/sqrt(Hz) to -1.0 dB at 30 fT/sqrt(Hz) (crossing near
  12 fT/sqrt(Hz)). Equal channel counts or coverage do not explain these differences (review).

## 8. Realistic adult OPM vs Neuromag comparison (G2, NEW) — `opmsquid.g2`, `scripts/g2_adult_comparison.py`

Configuration: `configs/g2_adult.toml`; assumptions A-G2-*, D-G2-COND in the register. Report:
`results/g2/G2_report.md` (generated from the JSON by `scripts/g2_report.py`).

Design
* Arrays: Neuromag T3 at the sample recording's measured head position (channel sets mag, grad and
  combined, with the full cross-type covariance); OPM arrays `opm_matched` (97 matched sites,
  coverage control), `opm204` (channel-budget control) and `opm_dense` (212 sites, full system),
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
  resampling targets instead gives CIs about 4-8x narrower),
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

Results (`results/g2/g2_summary.json`, v2 arrays, run at 81168f3: refined head surface and cell
clearance; medians over the 7,661 targets with parcel-bootstrap 95 % CIs; OPM 15 fT/sqrt(Hz);
3-layer BEM)
* Noise validation: the empty-room model (brochure intrinsic noise + fitted room field) gives
  115 fT (magnetometers) and 21.4 fT/cm (gradiometers) against 115 fT and 20.2 fT/cm measured; the
  room field explains 94 % of the magnetometer and 1.6 % of the gradiometer empty-room variance.
  Brain noise (task baseline minus empty room, 1-40 Hz): gradiometers 37.1 fT/cm (calibrated),
  magnetometers 262 fT measured vs 192 fT predicted (0.73x; the independent cortical background
  under-predicts magnetometer noise, likely distant and non-cortical sources). Median brain-noise
  RMS is ~500-525 fT at the OPM sites and 192 fT at the SQUID magnetometers.
* Intrinsic noise only (no brain noise): Neuromag combined beats every OPM array (dense 0.76x
  [0.74-0.79], matched 0.53x; OPM higher for <= 5 % of targets); OPMs beat the gradiometers alone
  (2.3-3.3x). The ratio scales as 1/(OPM noise): break-even at 11.5 (dense) and 7.9 (matched)
  fT/sqrt(Hz) against Neuromag combined.
* With brain noise (intrinsic+brain): matched 97-site OPM 1.00x [0.98-1.02] Neuromag combined
  (higher for 52 % of targets and 43 % of parcels: a tie), OPM 204 1.12x [1.10-1.16], dense
  212-site OPM 1.13x [1.10-1.16] (higher for 98 % of targets and all parcels); vs gradiometers
  alone 1.12-1.32x, vs magnetometers alone 1.03-1.18x. Adding the room field changes little (1.00,
  1.13, 1.14x). After the external projection (rank n - 8): 0.90x [0.81-0.95], 1.07x, 1.07x
  [1.01-1.13] (dense higher for 62 % of targets); relative to the unprojected room-field condition
  the projection costs the matched array 10 % of its detectability, the dense array 5 % and
  Neuromag 0.5 %.
* Depth (dense vs combined, intrinsic+brain): 1.66x at 10-15 mm, 1.39x at 15-20 mm, 1.23x at
  20-25 mm, 1.13x at 25-30 mm, 1.07x at 30-35 mm and 1.02-1.05x below 35 mm; the matched array is
  1.20x at 10-15 mm and falls to 0.91-0.95x below 40 mm. By distance from the inner skull (dense):
  1.44x at 4-5 mm, 1.09x beyond 10 mm. By lobe (dense): frontal 1.23x, parietal 1.14x, temporal
  1.14x, occipital 1.09x, insula 1.07x, cingulate 1.03x. Without the 558 medial-wall targets
  (7.3 %): dense 1.14x, matched 1.01x. (v1's rise of the dense ratio again below 45 mm, to 1.25x,
  was produced by cells too close to the coarse head surface.)
* Metric dependence: best single channel (peak-channel SNR) dense OPM 0.89x Neuromag (its best
  channel is usually a gradiometer), mean-power SNR 1.06x; detectability 1.13x.
* Estimated covariance (plug-in, Ledoit-Wolf): 10 s of data cost Neuromag combined 15 % and the
  dense OPM 9 % of the oracle detectability (60 s: 3 % and 2 %), so the dense/combined ratio is
  1.21x (10 s) and 1.14x (60 s); matched 1.12x and 1.03x.
* Extended sources: the ratio barely depends on patch size (dense 1.14, 1.14, 1.13x for 5, 10,
  20 mm radius; matched 1.01x), although cancellation reduces the net moment to 0.78, 0.53 and
  0.33 of the scalar moment.
* Sensitivity (dense | matched vs combined, intrinsic+brain): OPM noise 7 -> 30 fT/sqrt(Hz):
  1.28 -> 1.01x | 1.07 -> 0.93x; correlated background (lambda 5, 10 mm): 1.14, 1.16x | 1.01,
  1.02x; background calibrated on magnetometers: 1.13x | 1.00x; head position (+/-5 mm, +/-5 deg,
  well fitted): 1.09-1.16x | 0.98-1.02x; OPM scalp gap 3 and 6 mm (same sites moved out): 1.08 and
  1.04x | 0.98 and 0.96x. Joint OPM noise x scalp gap: dense 0.93x (30 fT/sqrt(Hz), 6 mm) to 1.13x
  (15 fT/sqrt(Hz), 0 mm), matched 0.83-1.00x. Frequency bands (brain scale and room field
  recalibrated per band): dense 1.13x (1-10 Hz), 1.13x (8-30 Hz), 1.11x (30-80 Hz; 1.08x with a
  100-Hz first-order OPM response); matched 0.99-1.01x.
* Convergence: background grid vs every usable vertex changes the median log2 ratios by <= 0.007;
  the 1-layer BEM refined from 5,120 to 20,480 triangles by <= 0.0023; Neuromag 4-point vs
  accurate integration 0.6 % in the gains; OPM 10-mm cell vs point 0.03 % median (95th percentile
  0.2 %; at most 1.0 % in the matched and 2.8 % in the dense array) in the peak-channel field of
  each of the 1,000 convergence targets; oct-6 vs random full-resolution targets <= 0.036;
  whitening tolerance none. Head model: on the convergence subset the dense/combined ratio is
  1.129x (primary), 1.147x with the v1 5,120-triangle head surface and 1.127x with a 1-layer BEM
  (matched 1.004, 1.008, 1.013x); gains differ by 11-12 % between 3 and 1 layers but the ratios
  hardly at all, and the 3- and 1-layer models agree to 0.2 %. In v1 (cells up to 2.4 mm inside the
  coarse head surface) the same comparison gave 1.21x vs 1.13x: the difference was numerical, not
  the skull and scalp.
* Bridge to G1A (intrinsic noise, peak-channel SNR): the realistic OPM/Neuromag magnetometer
  peak-field ratio follows the sphere with the real standoffs (7 mm OPM, 29.8 mm median SQUID);
  the equal-SNR depth at eta = 3 is 29.6 mm (matched) and 32.1 mm (dense) vs 29.1 mm (sphere, real
  standoffs) and 27.7 mm (Jas). OPMs are ahead at every depth for eta <= 2.0 (matched) or 2.25
  (dense) and behind at every depth for eta >= 5.0 (matched) or 5.5 (dense). The idealised
  benchmark's large OPM advantage assumes sensor-noise-limited SNR; with modelled brain noise,
  which on-scalp sensors also see more strongly, it shrinks to about 1.1x overall and vanishes
  for deep sources.

## 9. Epilepsy relevance, adult (G4, NEW) — `opmsquid.ied`, `opmsquid.detection`, `opmsquid.localization`, `scripts/g4_*.py`

Configuration: `configs/g4_epilepsy.toml`; assumptions A-G4-* in the register. The pediatric part
is section 11; head motion and OPM slippage (a bounded secondary extension) are section 12.

Detection design (`scripts/g4_epilepsy_adult.py`)
* Recordings: the G2 noise model in the time domain (A-G4-TS), one realization per 30-s segment
  shared by Neuromag, the matched 97-site OPM and the dense 212-site OPM array (OPM 15
  fT/sqrt(Hz)); 1-40 Hz, 150 Hz after decimation. Identical events (source, strength, morphology,
  time) in every array: 72 locations stratified by depth (10-20, 20-30, 30-45, 45-70 mm) and
  orientation (frontal-heavy: 24 frontal, 14 temporal, 14 parietal, 12 cingulate, 6 insular and 2
  occipital locations), focal dipoles at 10-320 nAm with three spike-wave morphologies, and 10-mm patches;
  1,728 events, one every 2 s.
* Detectors, per array and Neuromag channel set: a known-source/onset oracle (per-trial
  false-positive probability 0.001) and a practical scanner that knows neither time nor source
  (three waveform templates x 716 cortical candidates distinct from the true sources). Whitener from
  10 min of baseline null data; thresholds for 1 and 0.2 false events per minute from 20 min of
  independent calibration null data, frozen; held-out null (20 min) gives 0.5-1.1 and 0-0.2 false
  events per minute. A hit is the scan statistic above threshold within +/-50 ms of the true spike
  peak.
* Statistics: events share locations (18 per depth band, each with 6 strengths x 3 morphologies),
  so the location is the unit. Paired comparisons count, per location, the events detected by only
  one system and test the per-location differences with an exact sign-flip test
  (`detection.sign_flip_p`); the event-level McNemar test is kept for reference only. S50
  intervals come from a bootstrap over locations; the Wilson bands in the figure are event-level
  and descriptive.

Detection results (v2 arrays, run at 81168f3; focal, three morphologies pooled; strength for 50 %
detection, S50, with a bootstrap over locations; practical detector at 1 false event per minute; the
72 locations never lie on the medial wall; p-values uncorrected)
* S50, Neuromag combined vs dense OPM vs matched OPM: 53 vs 35 vs 41 nAm at 10-20 mm, 89 vs 82 vs
  88 nAm at 20-30 mm, 132 vs 127 vs 137 nAm at 30-45 mm, 290 vs 271 nAm vs not reached at 45-70 mm
  (95 % intervals 44-66, 28-50 and 32-59 nAm at 10-20 mm, which overlap; the deepest band's
  upper limits lie beyond the tested 320 nAm). The oracle
  needs about half the strength (25 vs 16 vs 20 nAm at 10-20 mm): searching over time and
  sources costs roughly a factor 2. The intervals of separate systems overlap even where the paired
  comparison is clear; the paired S50 ratio (same location resamples for both systems) is the
  comparison: Neuromag / dense 1.51 [1.25-1.66] at 10-20 mm, 1.08 [1.00-1.20], 1.04 [0.99-1.12] and
  1.07 [1.00-1.16] in the deeper bands.
* Paired, by location (locations with more events detected only by OPM / only by Neuromag
  combined; exact sign-flip p): dense OPM 16/0, 6/2, 6/1 and 6/1 across the four depth bands
  (p < 0.001, 0.23, 0.28, 0.13); matched OPM 9/3, 5/4, 4/5 and 0/7 (p = 0.04, 0.85, 0.62, 0.016:
  worse in the deepest band). With the oracle: dense 16/0, 13/0, 7/1, 5/5 (p < 0.001, < 0.001,
  0.06, 1); matched 10/5, 4/6, 2/5, 2/6 (p >= 0.22).
* Run-to-run variability: an earlier v2 run whose dense array differed by one site (so every noise
  draw differed) gave dense 13/0, 8/1, 8/2 and 5/2 (p = 0.0002, 0.03, 0.08, 0.8) and matched 9/2,
  6/2, 4/5 and 1/8 (p = 0.14, 0.23, 1, 0.02). The dense 10-20 mm advantage is stable; the 20-30 mm
  location-level result is not (p 0.03 there, 0.23 here) and is not claimed; v1's dense advantage
  in every band (14/0, 10/2, 12/2, 13/0) is superseded.
* Sensitivity at 1 false event per minute: superficial (10-30 mm) 40-nAm spikes 0.17 (Neuromag),
  0.28 (matched OPM), 0.31 (dense OPM); deep (30-70 mm) 160-nAm spikes 0.35, 0.31, 0.38 (frozen
  1-per-minute thresholds).
* Consistency with G2 (`scripts/study_g4_vs_g2.py`): at the same 18 locations per band, the G2 detectability ratio dense /
  combined (intrinsic+brain+env) has medians 1.46, 1.15, 1.07 and 1.03, and the oracle S50 ratio is
  1.62, 1.36, 1.10 and 0.98: the same ordering and direction; S50 pools the detection curves of 18
  locations over a factor-2 strength grid, so it need not equal the median ratio.
* Detection agrees with G2: the full OPM system detects superficial spikes at lower strength
  (about a third lower at 10-20 mm) and the advantage fades with depth: by location it is
  established at 10-20 mm for the practical detector (and at 20-30 mm only for the oracle). A
  matched-site OPM array is ahead by location only at 10-20 mm (9/3, p = 0.04 uncorrected; 0.14 in
  the earlier run: not robust) and behind at 45-70 mm (0/7, p = 0.016 uncorrected; the same
  direction in the earlier run).

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
  focal events: 66 % matched OPM, 57 % dense OPM, 38 % Neuromag). GOF, chi2/dof and the 95 %
  confidence volume are descriptive.
* Paired OPM-minus-Neuromag comparisons on identical events (one per location): median error
  difference with a bootstrap CI and Wilcoxon signed-rank p; exact McNemar p for joint detection +
  localization within 10 mm; 16 comparisons per OPM array, uncorrected.

Localization results (v2 arrays, run at 81168f3, `results/g4/g4_localization_summary.json`;
Neuromag combined, matched OPM, dense OPM)
* Detection in this subset: 320-nAm focal 0.92, 0.92, 0.92; 320-nAm patches 0.83 for every array;
  80-nAm patches 0.08-0.17 (too weak: their "localization" is that of noise).
* ECD, detected events: 320-nAm focal 4.3, 4.8, 4.2 mm on the MRI (2.8, 4.6, 2.5 mm in the sensor
  frame; the coregistration error alone displaces the source by a median 2.9 mm); 320-nAm patches
  7.2, 8.5, 6.7 mm; 80-nAm focal 5.1, 6.7, 3.7 mm. chi2/dof 0.9-1.1: the head-model and
  coregistration mismatch is small against the noise at these SNRs.
* dSPM peak, detected events: 320-nAm focal 13.6, 11.1, 11.2 mm; 320-nAm patches 18.5, 13.5,
  10.1 mm. Support recovery of the 320-nAm patches (share of the N strongest grid sources inside the
  patch; chance ~N/3,821): medians 0.1, 0.1, 0.3.
* Detected and localized within 10 mm (dSPM | ECD): 320-nAm focal 0.21 | 0.71, 0.38 | 0.62,
  0.33 | 0.75; 320-nAm patches 0.17 | 0.58, 0.38 | 0.58, 0.42 | 0.58.
* Paired (24 locations; 16 comparisons per OPM array, uncorrected p): for 320-nAm patches both OPM
  arrays localize better with dSPM (median error -4.8 mm matched, p = 0.14; -5.3 mm dense,
  p = 0.046) and more patches are both detected and localized within 10 mm (5/0, p = 0.06; 6/0,
  p = 0.03); none of these survives a correction over 16 comparisons. An earlier v2 run (one dense
  site fewer, hence different noise draws) gave -4.8 mm (p = 0.003) and -3.7 mm (p = 0.01): the
  direction and size are stable, the significance is not. No ECD difference and no focal-source
  difference was detected; differences are not excluded.
* With 24 locations the study bounds rather than resolves localization differences. In this model
  the dipole fit is limited by the coregistration error (about 3 mm) for every array, while the
  distributed estimate of extended sources tends to benefit from on-scalp sensors, matched or
  dense (about 5 mm; not established at this sample size).

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

G3B results (`results/g3b/g3b_summary.json`, `G3B_report.md`; computed at cb1b8a9 with all six
anatomies and redrawn at fd40dfe; the four anatomies of the independently reviewed pass at 8632ec9
reproduce to 1e-9, apart from 188 bootstrap bounds of secondary comparisons, which moved because
the added templates share the random-number stream (no headline interval changed). Dense OPM vs
Neuromag, intrinsic + brain noise, primary placement unless stated; dB of detectability;
area-weighted medians without the medial wall, parcel-bootstrap 95 % intervals. "Template" alone
means the 2-year template; the 18- and 12-month templates are named.)
* Heads: occipitofrontal circumference 586 mm (adult), 525 (school-age size, x0.895), 495 (2-year
  size, x0.844) and 495 mm (template); usable cortex 1,878, 1,470, 1,290 and 1,062 cm^2. Refitted
  OPM arrays: dense 212, 171, 155 and 151 sites (23.5-25.9 per 100 cm^2 of covered scalp, which
  shrinks from 901 to 584 cm^2), matched 96, 89, 90 and 83. Top contact raises the heads by 5.5,
  22, 22 and 28 mm; the median magnetometer-to-scalp gap is then 28.4, 34.1, 39.1 and 38.5 mm
  (centred 29.8, 41.3, 46.9, 47.3 mm). Every placement and the counterfactual helmets are feasible;
  the counterfactual factors are 0.895, 0.844 and 0.904 (template: under the adult's pose it sits
  right of the helmet's midline, centred left/right temporal gaps 48.9/33.3 mm, and its right side
  stops the helmet from shrinking to 0.844). Lateral centring moves the heads by -1.5 (adult), -0.5,
  -0.5 and -6.5 mm (template) along device x; about the laterally centred template the factor is
  0.864 (median gap 27.0 mm, adult 29.8 mm). The 18- and 12-month templates: circumference 491 and
  469 mm, usable cortex 975 and 895 cm^2, dense 157 and 145 sites (26.5 and 25.8 per 100 cm^2 of
  593 and 562 cm^2 covered scalp), matched 80 and 83; top contact raises them by 28 and 33.5 mm
  (median gap 38.9 and 40.4 mm; centred 46.0 and 48.9 mm). They too sit right of the midline
  (centred left/right temporal gaps 47.5/35.5 and 48.1/39.4 mm; lateral centring -5.5 and -2.5 mm),
  so their counterfactual factors are 0.893 and 0.851 (head-circumference ratios 0.838 and 0.801)
  and 0.863 and 0.841 about the laterally centred head.
* Link to G2: at the adult's measured position the dense/combined ratio is 1.130x as G2 reports it
  (unweighted, all targets), 1.143x without the medial wall and +1.25 dB area-weighted; at top
  contact the adult's D is +0.99 dB [+0.75, +1.20] (1.12x).
* D_child, D_adult, Delta (vs Neuromag combined): school-age size +1.49 vs +0.99 dB, Delta +0.47
  [+0.38, +0.58] dB; 2-year size +2.15 vs +0.99, Delta +1.13 [+0.99, +1.25]; 2-year template +1.85
  vs +0.99, Delta +0.73 [+0.44, +1.16] (66 parcels, Delta > 0 in 96 % of their area; the frontal
  poles are sparse). Against the gradiometers alone Delta is +0.69, +1.53 and +1.05 dB, against the
  magnetometers alone +0.44, +1.03 and +0.82 dB. With the external projection +0.90, +1.58 and
  +0.74 dB; for the matched-site OPM array +0.46, +1.22 and +0.90 dB (the gain is not a coverage
  effect of the dense array); for 5- and 10-mm patches +0.61/+0.42, +1.34/+1.14 and +0.58/+0.50 dB.
  Secondary metrics: peak-channel SNR (D_adult -1.19 dB, Neuromag's best channel ahead) Delta
  +0.77, +1.76, +1.65 dB; mean-power SNR (D_adult +0.55 dB) Delta +0.18, +0.58, +0.57 dB.
  18- and 12-month templates: D_child +1.90 and +2.06 dB, Delta +0.85 [+0.49, +1.27] and +1.03
  [+0.66, +1.44] dB (66 and 65 parcels, Delta > 0 in 96 and 94 % of their area); against the
  gradiometers +1.21 and +1.12, the magnetometers +0.67 and +1.01 dB; projected +0.91 and +1.17 dB;
  matched-site array +0.85 and +1.07 dB; 5- and 10-mm patches +1.07/+0.73 and +1.08/+0.94 dB;
  peak-channel SNR +1.36 and +1.62, mean-power SNR +0.77 and +0.92 dB. Across the three templates
  Delta grows as the head gets smaller (24, 18, 12 months: +0.73, +0.85, +1.03 dB), with
  overlapping intervals.
* By depth (template vs adult, combined, native depth strata): Delta +0.12 to +0.35 dB down to 30 mm
  (intervals include 0), +0.33 [+0.17, +0.60] at 30-40 mm, +0.37 at 40-50 mm, +1.12 at 50-60 mm and
  +2.23 at 60-90 mm (83 template vs 26 adult targets, 79 of them isthmus cingulate). Within strata the
  template's Delta is largest deep and for radial sources, but those strata hold little of its area
  (2.9 % of it is deeper than 50 mm); most of its pooled gain reflects its shallower cortex (42 % of its area at
  10-20 mm vs 21 % in the adult, area-weighted median depth 21.8 vs 26.2 mm): reweighted to the
  adult's depth mix, the target-level difference of the medians falls from +0.85 to +0.33 dB. For the
  scaled controls the homologous (vertex-wise) Delta is largest near the surface and falls with
  depth: 2-year size +2.25 [+1.67, +2.53] dB at an adult depth of 10-15 mm, +1.56 at 20-25 mm,
  +0.68 at 30-40 mm, +0.45 at 40-50 mm; school-age size +0.73, +0.70, +0.31 and +0.22 dB. By
  orientation: the scaled controls gain about equally for radial and tangential sources (+0.70 /
  +0.47 dB; +1.23 / +1.24 dB), the template far more for radial sources (0-30 deg +2.22 [+1.58,
  +2.82], 30-60 deg +1.00, 60-90 deg +0.53 dB; radial sources are 16 % of its targets). These are
  pooled over depth; at matched depth radial sources still gain more down to 40 mm (+0.64, +0.72,
  +0.81 vs tangential +0.21, +0.19, +0.09 dB at 0-15, 15-25, 25-40 mm) but not deeper (+0.48 vs
  +0.70 dB at 40-90 mm). The 18- and 12-month templates repeat this (`template_depth_checks`):
  within strata Delta is +0.21 to +0.83 dB at 10-50 mm (intervals exclude 0 from 25 mm down in both,
  and at 15-25 mm for the 18-month template) and +1.07 to +2.15 dB deeper; reweighted to the adult's
  depth mix their pooled difference falls from +0.91 to +0.54 and from +1.07 to +0.32 dB (36 and 46 %
  of their area at 10-20 mm; median depth 23.2 and 20.6 mm); at matched depth radial sources gain
  more than tangential ones down to 40 mm (18 months +1.06, +1.38, +1.07 vs +0.37, +0.38, +0.19 dB;
  12 months +1.36, +0.71, +0.91 vs +0.50, +0.30, +0.08 dB) and not deeper (+0.51 vs +0.62 and +0.67
  dB).
  Regionally the template's Delta is asymmetric (left inferior temporal, entorhinal, fusiform and
  pars orbitalis +2.3 to +3.0 dB; right orbitofrontal -0.8 to -0.9 dB): at top contact under the
  adult's measured pose its left frontal, temporal and parietal gaps are about 8 mm wider than the
  right ones (occipital 1.7 mm); the x-5mm variant, which roughly re-centres it, changes the median
  D by -0.1 dB.
* What drives Delta (each child placement against the adult at the same rule): with the child left
  centred (ear line where the adult's was) Delta is +1.15, +1.73 and +1.79 dB; at top contact +0.47,
  +1.13 and +0.73 dB; laterally centred, then top contact, +0.47, +1.12 and +0.68 dB; at true 18-mm
  contact +0.46, +1.16 and +0.70 dB; back contact +0.89, +1.30 and +1.39 dB. In the counterfactual
  helmet scaled with the head it is -0.17 [-0.21, -0.12], -0.30 [-0.35, -0.24] and +0.28 [+0.06,
  +0.53] dB, and about the laterally centred head -0.17, -0.30 and -0.37 [-0.54, +0.11] dB: the
  template's positive residual came from its lateral offset. For the 18- and 12-month templates:
  centred +1.93 and +2.49, top +0.85 and +1.03, laterally centred then top +0.69 and +0.99, 18-mm
  contact +0.83 and +0.95, back +1.63 and +1.68 dB; counterfactual +0.20 [-0.05, +0.58] and +0.08
  [-0.23, +0.22] dB, about the laterally centred head -0.29 [-0.45, +0.01] and -0.12 [-0.41, +0.10]
  dB. The relative gain of the head-adaptive array in children therefore comes from the fixed
  helmet's fit; with a helmet that fits as the adult's does, it vanishes or reverses slightly.
* Mechanism. Both systems' detectability rises in the smaller heads, by different routes:
  vertex-wise from the adult, the dense OPM gains +1.15 and +1.74 dB in the two scaled controls,
  Neuromag combined +0.65 and +0.54 dB (gradiometers +0.45 and +0.10, magnetometers +0.69 and
  +0.64 dB). The on-scalp OPM sees more signal from a cortex closer in absolute terms (median peak
  field of a 10-nAm dipole 211 fT in the adult, 242, 262 and 271 fT in the children; 255 and 293
  fT at 18 and 12 months) while its brain noise stays near 450-540 fT. In the fixed helmet the
  cortex is farther from the SQUIDs (median target-to-magnetometer distance at 15-20 mm depth 48
  mm adult, 50-57 mm children; OPM 26 mm in every head) and, with the background fixed per unit area, the smaller cortex carries
  less background power, so their brain noise falls (magnetometers 203 -> 150, 125 and 120 fT, 110
  and 98 fT at 18 and 12 months; gradiometers 41 -> 30, 23, 23, 21 and 18 fT/cm against 21 fT/cm
  intrinsic) more than their signal (magnetometer peak 70 -> 66, 54, 65, 63 and 67 fT). Absolute
  detectability of a 10-nAm dipole (median dB): dense OPM -1.60 (adult), -0.66, -0.19, +0.93
  (template), +1.08 (18 months) and +1.87 (12 months); Neuromag combined -2.66, -2.15, -2.27,
  -1.19, -1.02 and -0.62. In the counterfactual helmet the SQUID gains at
  least as much as the OPM: the counterfactual Delta depends on the comparator (gradiometers -0.32
  and -0.52 dB, magnetometers -0.12 and -0.21 dB for the scaled controls), which points to a
  SQUID-side component (signal against fixed intrinsic noise, largest for the gradiometers, whose
  brain noise is only about twice their intrinsic noise) rather than to the OPM's fixed standoff
  alone; the noise composition in the counterfactual helmet is not reported.
* Channel count: the children's dense arrays have fewer sites (171, 155, 151). The adult's dense
  array subsampled (farthest-point) to those counts has D +0.66, +0.52 and +0.49 dB (full array
  +0.99), so at an
  equal channel count Delta would be +0.72, +1.50 and +1.36 dB (18- and 12-month templates, 157 and
  145 sites: adult subsampled +0.54 and +0.41, Delta +1.37 and +1.49 dB): the smaller site count of
  the head-adaptive array works against it, and the headline Delta includes that loss.
* Placements: the source-blind variants (+-5 mm, pitch +-10 deg, roll +-5 deg, back contact) give
  D_child 1.45-2.52 (school-age size), 1.70-3.06 (2-year size), 1.75-2.51 (template), 1.81-2.49
  (18 months) and 1.99-2.80 dB (12 months), adult 0.96-1.41 dB; top contact (the primary) ranks in
  the middle of the family. Among the placements, true 18-mm contact gives the lowest D for the
  adult, the school-age size control and the three templates (0.86, 1.35, 1.70, 1.75 and 1.86 dB),
  but not for the 2-year size control (2.05 dB; pitch +10 deg gives 1.70). Regions:
  raising the head (centred -> top) lowers D in every lobe, most in the parietal and frontal
  lobes (template frontal +3.71 -> +2.02, parietal +3.02 -> +1.39 dB); back contact favours the
  SQUID at the occiput (template occipital +1.25 dB vs +2.08 at top) and disfavours it frontally
  (+3.89 vs +2.02 dB).
* Sensitivity (difference of the median D, child minus adult, with the same variant applied to both;
  a check on the medians, not the paired Delta estimator; vs combined): OPM noise 7-30 fT/sqrt(Hz)
  +0.38 to +0.55 (school-age size), +1.03 to +1.20 (2-year size), +0.64 to +0.87 dB (template);
  background variance x0.5 / x2 +0.51 / +0.48, +1.20 / +1.10, +0.89 / +0.80 dB (the variants scale
  the adult too; D_child itself moves by at most 0.06 dB, template 1.86 / 1.79 vs 1.85 dB); 1-layer
  BEM +0.43, +1.10, +0.78 dB. 18- and 12-month templates: OPM noise +0.65 to +1.01 and +0.81 to
  +1.13 dB, background x0.5 / x2 +0.89 / +0.91 and +1.10 / +1.01 dB, 1-layer BEM +0.79 and +1.00 dB.
  At 30 fT/sqrt(Hz) the adult's D is -0.04 dB (a tie) and the children's +0.34, +0.99, +0.61, +0.61
  and +0.77 dB.
* Usefulness (d >= 5) at 100 nAm, share of usable cortical area, both / OPM only / SQUID only /
  neither: adult 0.66 / 0.02 / 0.00 / 0.32, school-age size 0.69 / 0.04 / 0.00 / 0.28, 2-year size
  0.69 / 0.05 / 0.00 / 0.26, template 0.75 / 0.05 / 0.00 / 0.20, 18 months 0.76 / 0.05 / 0.00 / 0.19,
  12 months 0.79 / 0.05 / 0.00 / 0.16; at 50 nAm OPM only 0.07, 0.10, 0.14, 0.12, 0.12 and 0.13. The
  SQUID alone is usable on at most 0.15 % of the area; the cortex neither reaches is mostly deep and
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

Pediatric G4 results (`results/g4/G4_pediatric_report.md`, `g4_pediatric_comparison.json`; 2-year
template simulated at 1b2cb02, school-age control at 2c5735d, 2-year size control at c950959, 18-
and 12-month templates at cb1b8a9, comparison at fd40dfe; the code paths of the anatomy, the arrays
and the studies are the same at these commits; 18 locations per depth band in every anatomy;
p-values uncorrected, over 24 paired detection comparisons (3 comparators x 2 detectors x 4 bands)
and 16 localization comparisons per OPM array and anatomy). The children's Neuromag is at the
primary G3B placement (top contact), the adult's at its measured position (the G4 adult study as
frozen); top contact would raise the adult's head by only 5.5 mm, against the children's 22-33.5
mm. Order below: school-age size, 2-year size, 2-year template, 18 months, 12 months.
* Held-out null: 0.4-1.7 false events per minute at the 1-per-minute thresholds (adult 0.5-1.1):
  thresholds calibrated on 20 min and checked on 20 min differ by up to this much, so operating
  points are approximate, and they are lopsided between arrays in either direction (dense OPM vs
  Neuromag combined: 1.0 vs 0.55, 0.9 vs 0.4, 1.25 vs 0.65, 1.05 vs 0.95, 0.70 vs 1.25 per minute).
  At a matched held-out rate of 1 per minute the superficial advantage remains (sensitivity for
  40-nAm spikes at 10-30 mm, dense OPM vs Neuromag combined: 0.31 vs 0.18, 0.31 vs 0.13, 0.43 vs
  0.23, 0.50 vs 0.38, 0.46 vs 0.32; adult 0.31 vs 0.17).
* Strength for 50 % detection, practical detector, Neuromag combined vs dense OPM, 10-20 mm: adult
  53 vs 35, then 52 vs 37, 59 vs 36, 45 vs 30, 36 vs 27 and 40 vs 29 nAm; 45-70 mm adult 290 vs 271,
  then 275 vs 245, 273 vs 243, 243 vs 229, not reached vs 282 and 279 vs 259 nAm. The templates'
  point estimates are lower than the adult's at 10-30 mm, but the intervals overlap (e.g. Neuromag
  at 20-30 mm 71 [61-96], 62 [54-76] and 64 [52-91] vs 89 [67-122] nAm), and the sampled locations'
  median depths differ between anatomies by a few mm (10-20 mm: 13.6-17.8 mm, adult 16.9; 45-70 mm:
  48.0-54.3 mm, adult 51.7): no difference between the anatomies is claimed.
* Paired, dense OPM vs Neuromag combined (locations favouring OPM / Neuromag; paired strength ratio
  Neuromag / OPM): at 10-20 mm every anatomy favours the OPM: adult 16/0 (1.51 [1.25-1.66]), then
  16/0 (1.42 [1.20-1.57]), 16/0 (1.64 [1.40-1.85]), 16/2 (1.50 [1.26-1.68]), 16/0 (1.33
  [1.16-1.62]) and 14/1 (1.39 [1.19-1.68]); p <= 0.0005 in each. Deeper, a location-level
  difference appears only for the 2-year size control (30-45 mm 10/1, p = 0.009, ratio 1.18
  [1.06-1.31]; 45-70 mm 8/1, p = 0.03, 1.12 [1.03-1.24]) and the 18-month template at 45-70 mm (7/1,
  p = 0.047; Neuromag's S50 not reached, ratio interval 1.03-1.29); none survives a correction over
  24 comparisons, and elsewhere p >= 0.09 (ratios 1.00-1.12; for the school-age control at 45-70 mm
  the ratio's interval excludes 1, 1.12 [1.04-1.22], while p = 0.11). With the oracle every anatomy favours
  the OPM at 10-20 mm (12/1 to 17/1), and most do at 20-30 mm and in one or both deeper bands (e.g.
  18 months 9/2, 10/3, 12/1; 12 months 10/0, 8/1, 7/1). The matched-site array favours the OPM at
  10-20 mm with the practical detector in the 2-year size control (12/1, p = 0.002) and the 18- and
  12-month templates (12/0, p = 0.0005; 9/1, p = 0.014), not established in the school-age control
  (6/3) or the 2-year template (11/5, p = 0.05); the adult's deficit in its deepest band (0/7) is not
  seen in any child.
* So the detectability gains of G3B (Delta +0.47 to +1.13 dB) are not resolved by the practical
  detector with 18 locations per band beyond the superficial band, whose advantage is about the same
  in the adult and every smaller head (strength ratio 1.3-1.6); they appear mostly in the oracle's
  deeper bands.
* Localization (24 locations; Neuromag, matched, dense): ECD errors of detected events are similar
  in every anatomy (320-nAm focal: 4.1-7.8 mm). dSPM, all events, 320-nAm focal: the dense array is
  paired-closer than Neuromag in the three templates (-4.7, -5.9 and -4.7 mm; p = 0.035, 0.002 and
  0.005) and the school-age control (-1.1 mm, p = 0.049), not established in the 2-year size control
  (-3.6 mm, p = 0.06) or the adult (0.0 mm); 320-nAm patches: 2-year template -1.2 mm (p = 0.008),
  12 months -6.1 mm (p < 0.001), adult -5.3 mm (p = 0.046). ECD, 80-nAm sources, 18 months: -3.5 mm
  (focal, p = 0.003) and -4.0 mm (patches, p = 0.01). Detected and localized within 10 mm (dSPM),
  320-nAm focal spikes, Neuromag vs dense: adult 0.21 vs 0.33; children 0.25 vs 0.29, 0.12 vs 0.17,
  0.29 vs 0.54, 0.29 vs 0.58, 0.17 vs 0.62. None of these survives a correction over 16 comparisons
  per anatomy, and in the adult the significance of such differences varied between runs (section
  9); the consistent direction across the three templates (dSPM of strong focal events about 5 mm
  better with the dense OPM) is the more robust observation.
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
  out along their axes to the A-OPM-CLEAR clearances (a cap can lift, not sink). For each case,
  'known': the detectability with the displaced geometry (template and noise of the displaced
  geometry; the ideal limit of movement compensation or of a measured slip); 'mismatched': the
  output SNR (h^T s_k) / sqrt(h^T C_k h) of the matched filter h = C_k^+ s_0 built with the
  template s_0 of the reference geometry and the noise covariance C_k of the displaced data
  (`motion.mismatched_detectability`; a wrong-polarity template counts as -60 dB).
* B. In-band motion in a static residual field. A rigid array moving with the head reads
  n_i(t) . B(p_i(t)) at every integration point, with B(x) = B0 + G (x - x_ref) static in the room
  (`motion.readings`, exact for any rotation; `motion.jacobian`, its linearisation in the rotation
  vector about a pivot and the translation). With perfect calibration, rotation in B0 and
  translation in G change the readings by a uniform field in the head frame (removed exactly by a
  homogeneous-field projection, for any rotation), and rotation in G adds a symmetric traceless
  gradient (removed to first order by the 8-term projection of G2); `tests/test_motion.py` checks
  these identities, the Jacobian against finite differences and the basis against
  `environment.external_basis`. Calibration errors (A-MOT-CAL: sensitive-axis tilt 0, 1, 3 deg RMS
  with gain errors 0, 1, 3 % RMS, unknown to the analyst) leave residuals proportional to the field
  change. Head rotation about a pivot 60 mm below the head origin (A-MOT-PIVOT; the head origin as a
  sensitivity), independent per axis with in-band RMS theta; fields of unit strength (A-MOT-FIELD:
  1 nT uniform, or a 1 nT/m symmetric traceless gradient with unit Frobenius norm; 8 isotropic
  draws with their own calibration draws per case). The artefact covariance theta^2 P J J^T P^T
  after the correction P (none, the 3-term homogeneous or the 8-term projection, applied to signal
  and noise alike) is evaluated two ways (A-MOT-METRIC): outside the analyst's noise model (the
  matched filter of the static covariance applied to data that contain the artefact), and inside
  it (the oracle covariance: a rigid array's linear motion artefact spans at most three spatial
  patterns per field, which the optimal filter nulls; the bound for data-driven nulling or for
  regression on measured head motion). The rotation at which the median OPM detectability falls by
  1 or 3 dB, or to Neuromag's static detectability (D = 0), is interpolated on a logarithmic grid
  (0.0001-5 deg RMS); everything is linear in rotation x field, so for a field of B nT the
  thresholds divide by B. The SQUIDs are fixed in the room and have no such term.
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
computed at 5bfe1c0; medians over cortical targets of the adult and the 24- and 12-month templates;
dB of detectability, dense OPM or Neuromag combined, intrinsic + brain noise)
* A. Neuromag, head displaced in the fixed helmet, template of the reference position: 2 mm costs
  0.04-0.14 dB, 5 mm 0.28-0.56 dB, 10 mm down 1.65-1.72 dB (5-8 % of the cortex losing more than 3
  dB) and 10 mm sideways or forward/back 1.37-1.46 dB (templates only: at top contact the adult's
  head has no room); 5-deg rotations 0.22-0.84 dB, 10-deg rotations 0.91-2.31 dB (pitch and roll of
  the templates 1.90-2.31 dB, with 34-41 % of their cortex losing more than 3 dB). With the displaced
  geometry known (ideal movement compensation) the loss is at most 0.39 dB (10 mm away from the
  helmet top). 18 of the 81 displacements are infeasible at top contact (up or towards the wall).
  Dense OPM cap slipped, template of the reference geometry: 1 deg (median sensor shift 1.3-1.9 mm)
  costs 0.01-0.09 dB, 3 deg (3.8-5.7 mm) 0.11-0.53 dB; with the slip known at most 0.08 dB. Per mm of
  sensor-to-head displacement the uncompensated losses of the two systems are similar; the
  head-mounted array is unaffected by head displacement itself, which costs the fixed helmet up to
  1.7 dB at 10 mm and 2.3 dB at 10 deg without compensation.
* B. In-band rotation, artefact outside the noise model; thresholds in deg RMS per axis for a unit
  field (divide by the residual field in nT or nT/m), 1-dB loss, ranges over the three anatomies:
  no correction 0.021-0.023 deg (uniform) and 0.17-0.21 deg (gradient); homogeneous projection:
  the uniform term is removed exactly with perfect calibration, 0.43-0.45 deg with 1-deg/1-%
  calibration errors and 0.14-0.15 deg with 3 deg/3 %, while the gradient term is not removed
  (0.14-0.21 deg, as without correction); 8-term projection: uniform 0.38-0.41 and 0.13 deg,
  gradient 2.9-4.1 and 1.1-1.3 deg (1 deg/1 % and 3 deg/3 %). The OPM falls to Neuromag's static
  detectability (D = 0) at similar rotations (static D +0.51 to +1.86 dB after the projections). With
  the head origin as pivot the gradient thresholds change by less than 7 % after the homogeneous
  projection and rise by 10-33 % (or beyond 5 deg) after the 8-term projection. With the
  artefact part of the noise model (oracle), the loss stays below 0.15 dB up to 5 deg in every case.
  Artefact per channel for 1 deg RMS in the unit field: 14-16 pT without correction (uniform), 0.23-
  0.24 pT after the homogeneous projection with 1 deg/1 % errors, 22-29 fT for the gradient term
  after the 8-term projection. For example, with 1-deg/1-% calibration and the 8-term projection the
  dense array loses 1 dB at about 0.4 deg RMS of in-band head rotation in a 1-nT residual field, 0.04
  deg in 10 nT; without any correction at 0.02 deg in 1 nT.
* B'. Exact rigid motion (adult; 5-deg drift, 0.05-deg in-band jitter, 2 nT and 5 nT/m, 1 deg/1 %):
  the linear model predicts the exact in-band artefact to within 0.04 % (1,387, 365 and 21 fT RMS
  per channel after no, the homogeneous and the 8-term correction, against 89 fT of intrinsic OPM
  noise in the band); the drift moves each sensor's operating point by up to 274 pT (median 127 pT).
* Reading. A head-mounted array removes the geometry error that head motion causes in a fixed
  helmet, and converts motion into a field artefact whose size scales with the residual field and
  whose removal depends on calibration (or on modelling it from the data or from measured motion).
  These are bounds on two mechanisms in declared conditions, not an estimate of motion robustness
  in children.
