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
follows G3B.

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

## 10. Pediatric extension (G3) — `scripts/g3a_jas_size_benchmark.py`; G3B follows this baseline

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

### G3B: fixed adult helmet vs head-adaptive OPM (NEW)
Follows the adult baseline (owner decision 2026-09-30: the 2-year infant template and the adult
scaled to school-age and 2-year head size as size-only controls; no school-aged native anatomy).

## To be written

G3B; G4 pediatric.
