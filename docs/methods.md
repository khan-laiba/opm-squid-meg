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
a declared helmet-to-scalp gap (A-OPM-GAP), sensitive axis along the normal of the BEM head
surface (on the MRI scalp, A-BEM-CONFORM; normals from its own triangles) averaged within 15 mm
(A-OPM-AXIS; against a plane fitted to the dense MRI scalp within 15 mm of the site, a median
1.2 deg off, 95th percentile 4.9 deg matched and 7.4 deg dense, the largest at lower occipital
sites and next to the pinna, where no plane fits the scalp; v1-v3 took the adult's axes from its
outer skin's stored normals, a median 8-9 deg off that plane).
Clearance (A-OPM-CLEAR): every sensing centre is at least 6 mm (standoff - 1 mm) from the MRI
scalp and 4 mm from the BEM head surface, and the 10-mm cell is kept at least 1 mm outside the BEM
head surface at sample points: its 27 integration points and its surface, each face sampled at
11 x 11 points (1-mm spacing, corners and edges included; `opm.cell_volume_points`;
point-to-triangle distance over the triangles around the 10 nearest mesh vertices, the cell
oriented as the forward model builds it in the head frame). Every point of the cell surface lies
within 0.71 mm of a sample and the distance is 1-Lipschitz, so the rule itself guarantees only
0.29 mm for the whole cell (given exact distances at the samples) (with its centre outside and the head surface one closed surface far
larger than a cell, the cell is then outside it), not 1 mm. The exact cube-to-mesh distance of the
final cells (`opm.exact_cell_clearance`: the closest vertex-face and edge-edge pairs over every
triangle that could be closer, with an intersection test) is at least 1.056 mm (matched) and
1.006 mm (dense) in the adult (v4; recorded in `g2_summary.json`); over the nine G3B anatomies
every cell is at least 0.986 mm out (child A, matched array; 0.995 mm over the other eight, 18-month
template, dense array; `g3b_summary.json`, sensor distances).
Sites are moved outward along their axis in 0.5-mm steps by at most 5 mm, otherwise they are
infeasible and dropped. v2 checked only the
integration points, which reach +-3.87 mm of the +-5-mm faces: a check on 2026-10-01 found
cell corners up to 0.89 mm inside the head surface at 57 dense and 19 matched sites, so v3 applies
the rule to the whole cell: it moves most sites outward (median 1.0 mm from their v2 positions, the
matched sites by at most 2.5 mm; measured from the nominal standoff the moved sites sit 0.5-5 mm
out) and drops 7 dense and 2 matched sites, after which the greedy packing places a few dense sites
differently. In v3 the adult's BEM head surface lay about 1 mm outside its MRI scalp, so the rule
moved most of the adult's sites (median sensing-centre height 7.8 mm, vs 7.0 mm in the children,
whose head surfaces lie on their scalps); from v4 every anatomy's head surface lies on its MRI
scalp (A-BEM-CONFORM, section 3) and the rule moves only the sites the anatomy demands (adult: 11
of 98 matched and 26 of 208 dense sites; median height 6.99-7.01 mm in every anatomy and array). The rule
keeps the forward model's outer boundary out of the cell (see section 3, A-BEM-SKIN); the 27-point
cell average agrees with a 1,000-point one to 1e-3 (`tests/test_opm.py`).

Matched-site array: each Neuromag location is projected along its inward coil normal onto the
MRI scalp; sites below the brow plane (face/neck) or within 20 mm of a preauricular point are
dropped (A-OPM-COVER). Sample subject (v4): 98 of 102 sites (all that pass the coverage rules);
11 moved out by 0.5-5 mm (v3: 76 of 95; v2: 19); nearest-neighbour spacing 22.1-30.2 mm (median
26.5 mm); sensing centres 6.7-9.9 mm (median 7.0 mm) above the BEM head surface (Neuromag coils:
25-41 mm from the scalp, median 31).

Dense arrays (G2): farthest-point sampling (15-mm spacing) of the scalp in the same coverage
region, excluding scalp points more than 2 mm inside or 4 mm outside the smooth BEM head surface
(ear canals and pinna folds; the pinna) and within 20 mm of the lower edge of the MRI field of
view, then clearance as above and pruning to a 17-mm minimum centre spacing (A-OPM-PACK). The
field-of-view rule matters because the sample MRI is pitched ~36 deg relative to the head frame:
the flat cap where its head surface is cut lies above the brow plane at the back of the head.
The densest array found has 208 sites (`opm_dense`, v4; greedy, not proven maximal;
nearest-neighbour spacing 17.0-27.7 mm, median 19.2 mm; 26 sites moved out by up to 3.5 mm, v3: 164
by up to 5 mm, v2: 50; centres 6.5-11.0 mm, median 7.05 mm, above the BEM head surface); the
channel-budget control `opm204` takes 204 of them by farthest-point sampling, practically the full
system on this head (the channel-count control with Neuromag's three channels per site is the
triaxial array of section 8, A-OPM-TRIAX). The dense array reaches the lower occipital scalp
(lowest site, the sensing centre minus 7 mm along its axis: MRI z -105.5 mm vs -67.0 mm for the
matched array): v3's outer skin, 1-2 mm outside the scalp there, kept these cells from clearing it
within the 5-mm shift (v3's lowest dense site by this measure was at -72 mm; by its nearest MRI-scalp
point, the measure of section 8, at -75 mm). A 306-channel single-axis array appears infeasible under this rule (best spacing
found under the v2 rule 13.7-14.2 mm < 17 mm), though that is not proven. The composition and
placement of every array are pinned by `tests/test_g2.py`.

## 3. Forward models — `opmsquid.forward`, `opmsquid.anatomy`

MNE-Python 1.13.2 BEM forward models (single-compartment inner skull, 0.3 S/m, unless a paper
configuration prescribes otherwise; G2 uses the 3-layer 0.3/0.006/0.3 S/m model).

Outer conductor surface (A-BEM-CONFORM, v4). The OPM placement measures the standoff from the
MRI head surface, the forward model and the clearance rule use the BEM head surface, so the two
must agree in the same way on every head. The infant templates' BEM head surfaces are built on
their MRI head surfaces (each of their 2,562 vertices is a vertex of it). The sample subject's
stored outer skin, from its own segmentation, is not: over the OPM coverage region it lies a
median 0.80 mm outside the MRI scalp (90th percentile 1.5 mm, 99th 2.6 mm), and its stored vertex
normals are not those of its triangles (they lie within 5 deg of the radial direction from the
surface centroid). In v1-v3 this pushed the adult's OPM sensors outward under the whole-cell
clearance (164 of 205 dense sites moved; median height above the MRI scalp 7.78 mm, vs 7.00 mm in
every child) and tilted the adult's sensitive axes (A-OPM-AXIS): two adult-only differences in
G3B's D_adult. From v4 each vertex of the sample's head surface is moved at loading to the
nearest vertex of its MRI head surface (`anatomy.head_on_scalp`; same mesh and triangles; median
1.05 mm, 90th percentile 2.1 mm, at most 9.7 mm at the face and neck; refused if vertices merge,
a triangle flips or the outer skull leaves it; the outer skull stays at least 2.25 mm inside; no
self-intersection), its vertex normals following from its triangles. For the templates the rule
is the identity and the scaled controls inherit the conformed adult, so every anatomy's OPM
arrays sit at a median 6.99-7.01 mm above their MRI scalp. Every result that depends on the head
model or the OPM arrays was recomputed (`adult-baseline-v4`). Every 3-layer
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
`scripts/study_bem_skin_refinement.py`). On the v4 arrays (the head surface on the MRI scalp,
A-BEM-CONFORM: the cells sit at the nominal standoff, no integration point closer than 2.35 mm)
it gives a median 0.35-0.39 % at 3-4 mm (95th percentile 3.4-4.1 %) and 1.6-2.2 % at 2-2.5 mm
(13 points); the dense-array headline changes by -0.8 % (1.149x coarse, 1.141x refined; matched
1.008x both) and Neuromag gains by 0.10 %. More than 2 % separates the two surfaces on 7 dense and
4 matched channels (at most 4.4 % on the strong channels); the largest differences, for weak
fields (up to 0.65 of a source's array-RMS field at one channel), sit at the dense array's lower
occipital sites (MRI z -88 to -106 mm), where the conformed surface follows the occipito-cervical
scalp. The 1-layer model, which has no head surface, gives nearly the same headline (1.126x). An exact test on a
3-shell sphere shows MNE's BEM is accurate near a regular surface (95th percentile 0.45 % at 2 mm
with 5,120 triangles, 0.07 % with 20,480; `scripts/study_bem_sphere_accuracy.py`); the sample head
surface has larger triangles (median edge 7.2 mm, 95th percentile 13.6 mm) and a thin scalp in
places (the outer skull's vertices lie at least 2.25 mm inside the conformed surface, 1 % of them
within 2.9 mm). A second refinement is too large for a full BEM solution here. At the v3
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
  A-BEM-SKIN, which changes the Neuromag gains by 0.10 %; the JSON configuration field still reads
  5120; 0.33/0.0042/0.33 S/m, linear collocation); the sample recording's Vectorview geometry and head position; 4-point coil
  integration (study coils 9014/9024). The matched OPM array (98 sites) is the NEW column.
* Descriptors: depth to the nearest node of the 2,562-node BEM scalp; orientation to the normal of
  the nearest inner-skull BEM node, folded to 0-90 deg; the paper's 5 mm x 10 deg bins.
* Sources: 3783 dipoles sampled per bin to the paper's counts (the "dipole traces" cannot be
  rebuilt on new anatomy); 600 nAm PCHIP spike through the digitised Fig. 2(b) waveform. Patches
  grown from each dipole along mesh edges within +/-10 deg until the area first exceeds 20 mm^2
  (2921 of 3783 grow), uniform density with a 622-nAm median total (611-680 nAm; paper 612-678).
  The patch bins are not matched to the paper's counts (2,758 patches in its 7 x 8 grid vs 2,264);
  the 0-10 deg bins hold at most 7 patches, and a bin with fewer than 5 is left blank (below).
* Background: 28,494 random usable vertices (10 %), independent band-limited Gaussian moments,
  each peak-normalised to 10 nAm over a stationary 6 s (generated with 3 s of padding each side),
  one realization shared by all sources and arrays. Two levels (U-HU-bglevel): `as_specified`,
  and `fig6_calibrated`, which multiplies every background moment by 0.47: the paper's Fig. 6
  magnetometer baselines at channels 0631/0711/0741 divided by ours drawn and digitised
  identically, averaged over 20 independent background realizations because both the paper's
  baseline and one realization of ours are single draws (one realization alone gives 0.45-0.52,
  5-95 %; Appendix F of `docs/literature/hunold2016.md`). Other channel choices and statistics
  give 0.42 (gradiometers 0412/0413/0423) to 0.54 (all-channel median), which scales the
  calibrated bin means by 0.87-1.11. The scalar describes the effective background level at the
  sensors; whether the difference from the text's level comes from the normalisation, the
  filtering, the head position or the anatomy cannot be determined.
* SNR: per sensor type, the channel with the largest noise-free spike; background amplitude
  2 mean|hilbert| over the 1 s before onset (analytic signal of the whole trace, cropped; computed
  on the 1-s segment alone the amplitude is a median 0.2-0.4 % smaller, 2-3 % at the 5th
  percentile, and in the comparison below the mean ratios move by at most 0.005 and the 2.5
  classification by one patch bin: `hilbert_segment_only`);
  primary numerator the noise-free peak-to-peak spike (Appendix E: 0.87-1.10x the printed Fig. 6
  values depending on digitisation handling, against 0.61-0.73x for the literal noisy peak);
  variants peak, noisy peak and noisy peak-to-peak. Bin means as in the paper; unpaired per-bin
  tests (Student's t if both groups pass Shapiro-Wilk, else rank-sum). A patch bin with fewer than
  5 patches is left blank in the maps, the comparison and the tests (3 bins, 7 patches; the dipole
  bins follow the paper's counts).
* Comparison with the paper, per bin against the digitised colour classes of Figs 4(a) and 5(a)
  (midpoints). Calibrated, p2p: r = 0.96-0.98 (dipoles), 0.93-0.94 (patches); mean ratio
  0.74-0.85 (0.64-0.94 over the calibration range), bins with paper SNR >= 2.5 at 0.80-0.87x,
  weaker bins 0.70-0.80x; the 2.5 classification agrees in 79-89 % of bins; the GM - MM sign
  agrees in every patch bin and 53 of 54 dipole bins with a paper difference of at least one
  class. The noisy peak-to-peak is close to the paper on average (0.94-0.99x; strong bins
  0.85-0.94x, weak bins 0.99-1.08x; classification agreement 84-95 %), so the paper's values lie
  between the two numerators, closer to the noisy one; the primary numerator was fixed from
  Fig. 6 beforehand and is not changed to fit the maps. As specified, the maps keep their shape
  (same r) but reach only 0.35-0.40x the paper. v4 (the head surface on the MRI scalp,
  A-BEM-CONFORM) moved the sample's scalp nodes, from which depth is measured, about 1 mm inward:
  each depth bin now holds slightly deeper sources, which lowered the bin means by 5-7 %
  (v3: 0.78-0.91x).
* Fig. 6 spike check (Table 1 tangential examples; left posterior frontal cortex): the
  magnetometer spike agrees (superficial: ours 6.3 pT median, paper 6.4-7.3 pT) but our
  gradiometer spike is 1.5-1.7x the paper's (199 vs 119-134 pT/m); the gradiometer baseline is
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
  unaffected). The usable-vertex rule (A-BEM-DIST) also removes patch members within 4 mm of the
  inner skull: against the same patches over every valid vertex (median areas 2.86 and 7.54 cm^2)
  22 % (10 mm) and 31 % (16 mm) of the patches lose more than 5 % of their area and 8 % and 7 % more
  than 20 %; the pooled Eq. 1 SNR changes by a median 0 dB (5th percentile -1.3 and -1.5 dB; 7 % and
  9 % of patches by more than 1 dB; `patch_truncation` in `g1c_summary.json`). The same patch
  construction is used in G2, G3B and G4. G1C keeps the medial wall (aparc 'unknown': 5.9 % of the
  usable vertices, 558 of the 7,661 patch centroids, 97 of the 1,736 noise sources), as the paper's
  whole-cortex maps do, whereas G2 and G4 exclude it; without it the pooled focal median is -21.66
  instead of -22.13 dB and 58.9 instead of 56.0 % of the vertices lie in the paper's display range
  (`without_medial_wall`). With both changes, patches over every valid vertex (untruncated) at the
  7,103 centroids off the wall (v4, `variant_no_wall_untruncated`), the pooled patch medians are
  -26.42 (10 mm) and -21.58 dB (16 mm) instead of -27.22 and -22.25 dB, 60.6 and 60.2 % inside the
  display range instead of 56.5 and 60.8 %; the untruncated members include vertices within 4 mm of
  the inner skull, where the BEM is less accurate, so the variant bounds the two side-effects
  rather than replacing the primary rows.
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
  grid, 1.77 nAm when normalised to the ~4,000 sources of a 7-mm MNE grid (paper 1.6-1.9 nAm);
  focal SNR median -22.1 dB (pooled, all usable vertices; -22.6 dB at the oct-6 centroids; 5-95 %
  -35.2 to -16.7 dB); 56 % of vertices inside the paper's -29 to -19 dB display range. The deep
  medial regions are darkest (67 % below -29 dB in the cingulate, parahippocampal, entorhinal and
  medial orbitofrontal cortex; cingulate median -31 dB), as in the paper's MEG maps, while the
  medial occipital and paracentral cortex are bright; the insula (lateral but deep) is also low
  (-30 dB). Patch 16 mm minus 10 mm: median 5.4 dB (5.7 dB in the mesial temporal lobe) vs 8.3 dB
  for area scaling of our patches (areas 2.81 and 7.29 cm^2) and 10 dB quoted by the paper for the mesial temporal lobe
  (modality unclear, A15); cancellation within the larger folded patches explains the shortfall
  from area scaling. Skull 0.06 vs 0.006 S/m (each with its own calibration): median +0.45 dB
  (pooled), 95th percentile of the absolute change 3.4 dB, largest for near-radial sources.
  Maps use the paper's colour limits.
* Extension (NEW): the matched 98-site OPM array. The brain-noise sources are calibrated on the
  recorded minus empty-room variance (instrument noise is 6 % of the recorded variance for
  magnetometers and 35 % for gradiometers), then intrinsic noise is added explicitly (SQUID
  brochure values; OPM swept 7-30 fT/sqrt(Hz)), 0.5-100 Hz. Brain noise only, the OPM array ties
  both SQUID sensor types (Eq. 1 median +0.1 dB vs magnetometers, -0.08 dB vs gradiometers).
  With intrinsic noise, OPM vs gradiometers is +1.3 to +2.3 dB (the gradiometer noise floor),
  and OPM vs magnetometers +0.2 dB at 7 fT/sqrt(Hz) to -1.0 dB at 30 fT/sqrt(Hz) (crossing near
  12 fT/sqrt(Hz)). Equal channel counts or coverage do not explain these differences (checked).

## 8. Realistic adult OPM vs Neuromag comparison (G2, NEW) — `opmsquid.g2`, `scripts/g2_adult_comparison.py`

Configuration: `configs/g2_adult.toml`; assumptions A-G2-*, D-G2-COND in the register. Report:
`results/g2/G2_report.md` (generated from the JSON by `scripts/g2_report.py`).

Design
* Arrays: Neuromag T3 at the sample recording's measured head position (channel sets mag, grad and
  combined, with the full cross-type covariance); OPM arrays `opm_matched` (98 matched sites,
  coverage control), `opm204` (channel-budget control) and `opm_dense` (208 sites, full system),
  all single-axis, 10-mm cell, 7-mm standoff, no extra scalp gap in the primary run (section 2.2);
  and, from v4, a channel-count control `opm_triax`: the matched sites with three orthogonal axes
  each (294 channels; A-OPM-TRIAX).
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
  the estimate, evaluated under the true covariance; from v4 the estimation data of all arrays come
  from one common realization: the same background sources and room-field coefficients through
  every array, intrinsic noise per array). Comparisons: median over targets of
  log2(d_OPM / d_SQUID) with 95 % CIs from a bootstrap over 70 groups, the 68 Desikan-Killiany
  parcels and the two medial-wall ('unknown') labels (neighbouring targets are correlated;
  resampling targets instead gives CIs about 5-11x narrower),
  the share of targets and of parcels where OPM is higher, by depth, distance from the inner
  skull, orientation, lobe and on the cortex. Targets on the medial wall (FreeSurfer 'unknown')
  are kept in every median and reported separately.
* Sensitivity: OPM noise; correlated background; background calibrated on magnetometers; SQUID
  head position (measured; +/-5 mm along each device axis; +/-5 deg pitch; well fitted, 20 mm
  from the nearest magnetometer); OPM scalp gap (0, 3, 6 mm); a joint OPM noise x scalp gap grid;
  plug-in covariance; Neuromag sensor noise from the measured empty-room spectrum instead of the
  brochure values (v4, A-G2-SQUIDMEAS); frequency bands 1-10, 8-30 and 30-80 Hz with a 100-Hz OPM
  response (`scripts/g2_band_sensitivity.py`). Except the joint grid, one factor at a time: the analyses
  show the dependence, they do not bound it. The scalp-gap variants move the primary OPM sites
  outward along their axes (v2; v1 rebuilt the arrays, which added sensors at larger gaps).
* Convergence: background on the 7-mm grid vs every usable vertex; 3-layer vs 1-layer BEM and the
  1-layer BEM refined to 20,480 triangles (1,000 targets); Neuromag 4-point vs accurate coil
  integration; OPM point vs 10-mm cell; oct-6 vs random full-resolution targets; whitening
  tolerance.
* Bridge to G1A: the realistic OPM/magnetometer peak-field ratio vs depth and the equal-SNR depth
  d_eq(eta) (intrinsic noise only, peak-channel SNR) against the sphere benchmark. d_eq comes from
  the depth-binned medians of the ratio, not from per-source crossings.

Results (`results/g2/g2_summary.json`, v4: the head surface on the MRI scalp, A-BEM-CONFORM; run
at e53bea8; medians over the 7,661 targets with parcel-bootstrap 95 % CIs; OPM 15 fT/sqrt(Hz);
3-layer BEM)
* Noise validation: the empty-room model (brochure intrinsic noise + fitted room field) gives
  115 fT (magnetometers) and 21.4 fT/cm (gradiometers) against 115 fT and 20.2 fT/cm measured; the
  room field explains 94 % of the magnetometer and 1.6 % of the gradiometer empty-room variance.
  Brain noise (task baseline minus empty room, 1-40 Hz): gradiometers 37.1 fT/cm (calibrated),
  magnetometers 262 fT measured vs 192 fT predicted (0.73x; the independent cortical background
  under-predicts magnetometer noise, likely distant and non-cortical sources). Median brain-noise
  RMS is ~495-525 fT at the OPM sites and 192 fT at the SQUID magnetometers. The model matches one
  scalar of the measured brain noise: per channel the model/measured ratio spans 0.16-1.90
  (gradiometers, 5th-95th percentile, median 0.88) and 0.13-0.78 (magnetometers, median 0.55),
  part of which is measurement noise (one good gradiometer has a negative measured brain variance).
* Not modelled (limitations of every G2-G4 comparison): non-cortical physiological fields
  (cardiac, ocular, muscle), which on-scalp OPMs, being magnetometers, would see at about the
  amplitude of the SQUID magnetometers and which the 8-term projection removes only in their
  uniform-plus-gradient part; the 1/f rise of real OPM noise at low frequencies (white OPM noise is
  assumed; the 30-fT/sqrt(Hz) end of the sweep bounds a uniformly worse sensor, not a coloured
  one); the single-axis arrays have 208 channels against Neuromag's 306 (the triaxial control
  below matches the count under a declared tangential-noise assumption).
* Intrinsic noise only (no brain noise): Neuromag combined beats every OPM array (dense 0.74x
  [0.72-0.77], matched 0.53x; OPM higher for <= 4 % of targets); OPMs beat the gradiometers alone
  (2.3-3.2x). The ratio scales as 1/(OPM noise): break-even at 11.1 (dense) and 8.0 (matched)
  fT/sqrt(Hz) against Neuromag combined.
* With brain noise (intrinsic+brain): matched 98-site OPM 1.01x [0.99-1.03] Neuromag combined
  (higher for 56 % of targets and 47 % of parcels: a tie), OPM 204 1.14x [1.12-1.16] (with 208 dense
  sites this channel-budget control is practically the full system), dense 208-site OPM 1.14x
  [1.12-1.17] (higher for all but 4 of the 7,661 targets and for every parcel); vs gradiometers
  alone 1.13-1.33x, vs magnetometers alone 1.04-1.19x. Adding the room field changes little (1.01,
  1.14, 1.14x). After the external projection (rank n - 8; for Neuromag it acts on all 306 channels
  jointly, so "magnetometers" or "gradiometers" after the projection are subsets of the jointly
  projected data, not standalone sensor-type systems): 0.95x [0.90-0.98], 1.11x, 1.12x
  [1.08-1.15] (dense higher for 84 % of targets); relative to the unprojected room-field condition
  the projection costs the matched array 6 % of its detectability, the dense array 3 % and
  Neuromag 0.5 %. In v3 the projected dense ratio was 1.05x [0.99-1.10], a tie (next bullet).
* What the head-model change did (v3 -> v4; `scripts/study_head_surface_effect.py`,
  `results/g2/head_surface_effect.json`): the shallow ratios hardly moved (1.60x and 1.36x at
  10-20 mm in both), the deep ones rose (intrinsic+brain below 45 mm 1.05-1.07x, v3 0.98-1.01x;
  projected 0.92-1.00x, v3 0.69-0.76x), and the projection costs the dense array less (median
  retained signal norm 0.63 vs 0.56). The 1-layer BEM, which has no head surface, moved the same
  way (dense projected 1.106x vs 1.050x in v3), so the change comes from the arrays, not from the
  conductor's outer surface. Decomposition (the G2 machinery with point gains; it reproduces
  v3 exactly: dense 1.111x | 1.048x, matched 0.997x | 0.887x, intrinsic+brain | projected): the
  conformed surface in the forward model with the v3 arrays changes nothing (1.111x | 1.047x); the
  v4 sites with v3-style axes give 1.130x | 1.073x (matched 1.008x | 0.908x); the v4 axes, the
  conformed surface's own normals instead of the stored, nearly radial ones, bring 1.144x | 1.115x
  (matched 1.009x | 0.948x). Without the 12 dense sites below v3's lowest (nearest scalp point below
  MRI z -75 mm) the v4 dense array gives 1.132x | 1.092x: the lower occipital coverage adds about
  1 % | 2 %. The largest part of the projected-condition change is the axis definition (deepest band,
  45-70 mm: 0.73x in v3, 0.82x with the v4 sites, 0.98x with the v4 axes, 0.90x without the low
  sites): nearly radial axes made the external-field patterns more alike those of deep sources, so
  the projection removed more of their signal.
* Depth (dense vs combined, intrinsic+brain): 1.60x at 10-15 mm, 1.36x at 15-20 mm, 1.22x at
  20-25 mm, 1.13x at 25-30 mm, 1.08x at 30-35 mm and 1.05-1.07x below 35 mm; the matched array is
  1.20x at 10-15 mm and 0.94-0.98x below 40 mm. By distance from the inner skull (dense): 1.42x at
  4-5 mm, 1.10x beyond 10 mm. By lobe (dense): frontal 1.22x, parietal 1.14x, temporal 1.16x,
  occipital 1.12x, insula 1.09x, cingulate 1.05x. Without the 558 medial-wall targets (7.3 %): dense
  1.15x, matched 1.02x. (v1's rise of the dense ratio again below 45 mm, to 1.25x, was produced by
  cells too close to the coarse head surface.)
* Metric dependence: best single channel (peak-channel SNR) dense OPM 0.90x Neuromag (its best
  channel is usually a gradiometer; OPM higher for 28 % of targets), mean-power SNR 1.04x;
  detectability 1.14x.
* Estimated covariance (plug-in, Ledoit-Wolf, from one common noise realization): 10 s of data cost
  Neuromag combined 15 % and the dense OPM 9 % of the oracle detectability (60 s: 3 % and 2 %), so
  the dense/combined ratio is 1.23x (10 s) and 1.16x (60 s); matched 1.12x and 1.03x. The spread over
  realizations is not reported.
* Channel-count control (triaxial OPM at the 98 matched sites, 294 channels, A-OPM-TRIAX; intrinsic
  + brain | projected): with equal noise on every axis 1.12x [1.10-1.15] | 1.17x [1.15-1.19] Neuromag
  combined, 1.12x | 1.23x the matched single-axis array; with the tangential axes at twice the noise
  1.05x [1.03-1.08] | 1.09x [1.08-1.11]. At the same sites the two added axes are worth 5-23 %; at
  Neuromag's channel count the OPM array is ahead under these noise assumptions.
* Extended sources: the ratio barely depends on patch size (dense 1.14, 1.15, 1.13x for 5, 10,
  20 mm radius; matched 1.01x), although cancellation reduces the net moment to 0.78, 0.53 and
  0.33 of the scalar moment.
* Sensitivity (dense | matched vs combined, intrinsic+brain): OPM noise 7 -> 30 fT/sqrt(Hz):
  1.30 -> 1.01x | 1.08 -> 0.93x; correlated background (lambda 5, 10 mm): 1.15, 1.17x | 1.01,
  1.02x; background calibrated on magnetometers: 1.15x | 1.01x; head position (+/-5 mm, +/-5 deg,
  well fitted): 1.10-1.18x | 0.99-1.02x; OPM scalp gap 3 and 6 mm (same sites moved out): 1.09 and
  1.04x | 0.99 and 0.96x; Neuromag noise from the measured empty-room spectrum (median in-band
  25.9 fT and 20.0 fT/cm vs the brochure's 20.7 fT and 21.3 fT/cm): 1.17x | 1.02x (projected 1.13x
  | 0.96x). Joint OPM noise x scalp gap: dense 0.92x (30 fT/sqrt(Hz), 6 mm) to 1.14x
  (15 fT/sqrt(Hz), 0 mm), matched 0.84-1.01x. Frequency bands (brain scale and room field
  recalibrated per band): dense 1.15x (1-10 Hz), 1.14x (8-30 Hz), 1.11x (30-80 Hz; 1.07x with a
  100-Hz first-order OPM response); matched 0.99-1.01x.
* Convergence: background grid vs every usable vertex changes the median log2 ratios by <= 0.005;
  the 1-layer BEM refined from 5,120 to 20,480 triangles by <= 0.0015; Neuromag 4-point vs
  accurate integration 0.6 % in the gains; OPM 10-mm cell vs point 0.01-0.02 % median (95th
  percentile 0.1 %; at most 0.5 % in either array) in the peak-channel field of each of the 1,000 convergence
  targets; oct-6 vs random full-resolution targets <= 0.029; whitening tolerance none. Head model:
  on the convergence subset the dense/combined ratio is 1.141x (primary), 1.149x with the
  5,120-triangle head surface and 1.126x with a 1-layer BEM (matched 1.008, 1.008, 1.016x); gains
  differ by 11-12 % between 3 and 1 layers but the ratios hardly at all, and the 3- and 1-layer
  models agree to 1.3 % (0.3 % in v3: the conformed surface is rougher; section 3). In v1 (cells up
  to 2.4 mm inside the coarse head surface) the same comparison gave 1.21x vs 1.13x: the difference
  was numerical, not the skull and scalp.
* Bridge to G1A (intrinsic noise, peak-channel SNR): the realistic OPM/Neuromag magnetometer
  peak-field ratio follows the sphere with the real standoffs (7.0 mm OPM, 29.8 mm median SQUID;
  v3: 7.8 mm OPM); the equal-SNR depth at eta = 3 is 29.6 mm (matched) and 31.8 mm (dense) vs
  29.0 mm (sphere, real standoffs) and 27.7 mm (Jas). OPMs are ahead at every depth for eta <= 2.25
  (both arrays) and behind at every depth for eta >= 5.0 (matched) or 5.5 (dense). The idealised
  benchmark's large OPM advantage assumes sensor-noise-limited SNR; with modelled brain noise,
  which on-scalp sensors also see more strongly, it shrinks to about 1.14x overall and to 1.05-1.07x
  for deep sources.

## 9. Epilepsy relevance, adult (G4, NEW) — `opmsquid.ied`, `opmsquid.detection`, `opmsquid.localization`, `scripts/g4_*.py`

Configuration: `configs/g4_epilepsy.toml`; assumptions A-G4-* in the register. The pediatric part
is section 11; head motion and OPM slippage (a bounded secondary extension) are section 12.

Detection design (`scripts/g4_epilepsy_adult.py`)
* Recordings: the G2 noise model in the time domain (A-G4-TS), one realization per 30-s segment
  shared by Neuromag, the matched 98-site OPM and the dense 208-site OPM array (OPM 15
  fT/sqrt(Hz)); 1-40 Hz, 150 Hz after decimation. Identical events (source, strength, morphology,
  time) in every array: 72 locations stratified by depth (10-20, 20-30, 30-45, 45-70 mm) and
  orientation (frontal-heavy: 24 frontal, 14 temporal, 14 parietal, 12 cingulate, 6 insular and 2
  occipital locations), focal dipoles at 10-320 nAm with three spike-wave morphologies, and 10-mm patches;
  1,728 events, one every 2 s.
* Detectors, per array and Neuromag channel set: a known-source/onset oracle (per-trial
  false-positive probability 0.001) and a practical scanner that knows neither time nor source
  (three waveform templates x 716 cortical candidates distinct from the true sources). The
  scanner's templates are the three simulated morphologies, so its absolute sensitivity and
  false-event rates are optimistic; the paired comparison between arrays is less affected. The
  candidates' topographies come from the truth forward model (3-layer BEM, exact coregistration), a
  convenience shared by every array that makes the scanner optimistic in the same way. Whitener from
  10 min of baseline null data; thresholds for 1 and 0.2 false events per minute from 20 min of
  independent calibration null data, frozen; held-out null (20 min) gives 0.65-1.3 and 0.1-0.6 false
  events per minute (v4; v3 0.65-1.15 and 0.1-0.35); each 1-per-minute threshold rests on about 20 calibration events, whose Poisson
  95 % interval (0.6-1.5 per minute) matches that held-out spread, and the 0.2-per-minute one on
  about 4. A hit is an emitted event (a local maximum of the scan statistic, events at
  least the 0.25-s refractory period apart, as counted on the null data) above threshold within
  +/-50 ms of the true spike peak (v3, after a check; v2 took the statistic's maximum within the
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
  the first v3 summaries dropped resamples in which neither system reached 50 %, found by a check).
  The Wilson bands in the figure are event-level and descriptive.

Detection results (v4 arrays, simulated at ed852b2; focal, three morphologies pooled; strength for 50 %
detection, S50, with a bootstrap over locations; practical detector at 1 false event per minute; the 50 %
level and the 1-per-minute operating point are operational study choices, not clinical standards; the
72 locations never lie on the medial wall; p-values uncorrected)
* S50, Neuromag combined vs dense OPM vs matched OPM: 47 vs 34 vs 40 nAm at 10-20 mm, 86 vs 72 vs
  84 nAm at 20-30 mm, 139 vs 128 vs 143 nAm at 30-45 mm, 302 vs 254 vs 320 nAm at 45-70 mm
  (95 % intervals 36-63, 28-51 and 32-58 nAm at 10-20 mm, which overlap; in the deepest band
  Neuromag's and the matched array's upper limits lie beyond the tested 320 nAm). The oracle needs
  about half the strength (25 vs 16 vs 22 nAm at 10-20 mm): searching over time and sources costs
  roughly a factor 2. The intervals of separate systems overlap even where the paired comparison is
  clear; the paired S50 ratio (same location resamples for both systems) is the comparison:
  Neuromag / dense 1.36 [1.10-1.58] at 10-20 mm, 1.20 [1.06-1.34] at 20-30 mm, 1.09 [1.00-1.22] at
  30-45 mm and 1.19 [1.00-open] at 45-70 mm (Neuromag does not reach 50 % within the tested
  strengths in a third of that band's resamples).
* Paired, by location (locations with more events detected only by OPM / only by Neuromag
  combined; exact sign-flip p): dense OPM 11/1, 9/1, 5/1 and 9/0 across the four depth bands
  (p = 0.004, 0.014, 0.19, 0.004); matched OPM 7/4, 4/5, 2/5 and 3/6 (p = 0.40, 1, 0.77, 0.40). With
  the oracle: dense 15/0, 11/1, 8/2, 9/5 (p < 0.001, 0.008, 0.05, 0.14); matched 9/3, 5/5, 5/3, 4/6
  (p >= 0.09).
* Scoring and run-to-run variability: an array's channel count changes the random stream, so every
  array change also redraws the noise. v3 (whole-cell clearance, the adult's sensors pushed outward
  by its outer skin) gave dense 8/2, 5/1, 4/0 and 0/2 (p = 0.037, 0.19, 0.12, 0.50; S50 ratio 1.29
  [1.03-1.51] at 10-20 mm) and matched 8/4, 4/6, 1/3 and 1/9 (p = 0.52, 0.81, 0.62, 0.016); v2, which
  scored injected spikes by the statistic's maximum near the true peak, gave dense 16/0, 6/2, 6/1 and
  6/1 (p < 0.001, 0.23, 0.28, 0.12; S50 ratio 1.51 [1.25-1.66]) and matched 9/3, 5/4, 4/5 and 0/7;
  an earlier v2 run gave dense 13/0, 8/1, 8/2 and 5/2 and matched 9/2, 6/2, 4/5 and 1/8. The dense
  10-20 mm advantage holds in every run. v4 also favours the dense array at 20-30 mm and 45-70 mm.
  At the G4 locations the equal-standoff change raised the G2 ratio in every band, most in the
  deepest (by 1.5, 2.2, 3.2 and 4.9 % from 10-20 to 45-70 mm), yet the 30-45 mm band shows no
  detection difference (5/1, p = 0.19); of the earlier runs only the early v2 run found a deeper
  difference (20-30 mm, p = 0.03; otherwise 20-30 mm p = 0.23 and 0.19, 45-70 mm p = 0.12 and 0.50 in
  v2 and v3), so with 16 such tests per run the deeper results are reported, not established. The matched array's deficit at 45-70 mm (v3 and v2: 1/9 and 0/7, p = 0.016) is not
  seen in v4 (3/6, p = 0.40). v1's dense advantage in every band (14/0, 10/2, 12/2, 13/0) is
  superseded.
* Sensitivity at 1 false event per minute: superficial (10-30 mm) 40-nAm spikes 0.23 (Neuromag),
  0.29 (matched OPM), 0.36 (dense OPM); deep (30-70 mm) 160-nAm spikes 0.34, 0.31, 0.41 (frozen
  1-per-minute thresholds).
* Matched operating points (`scripts/study_g4_matched_rate.py`; the stored simulations
  re-evaluated with every detector's threshold set on the held-out null to 1 false event per
  minute, in-sample for those data; the re-summarised frozen-threshold values reproduce the
  committed ones exactly): dense OPM 12/1, 12/0, 5/1 and 9/0 (p = 0.002, 0.0005, 0.19, 0.004; S50
  ratio 1.40 [1.12-1.59] at 10-20 mm); matched OPM 8/4, 8/3, 2/5 and 5/6 (p = 0.32, 0.23, 0.77,
  0.79). The bands with p < 0.05 are the same as at the frozen thresholds.
* Consistency with G2 (`scripts/study_g4_vs_g2.py`): at the same 18 locations per band, the G2
  detectability ratio dense / combined (intrinsic+brain+env) has medians 1.45, 1.16, 1.09 and 1.05,
  and the oracle S50 ratio is 1.55, 1.29, 1.09 and 1.20: both favour the dense array most at
  10-20 mm and fall with depth to 30-45 mm; in the deepest band the oracle ratio rises again (9/5
  locations, p = 0.14) while G2's is 1.05. S50 pools the detection curves of 18 locations over a
  factor-2 strength grid, so it need not equal the median ratio.
* Detection agrees with G2: the full OPM system detects superficial spikes at lower strength
  (about a quarter lower at 10-20 mm) and the advantage fades with depth, with a deeper advantage
  in v4 that agrees in direction with G2's deep ratios but is not established across runs; a matched-site OPM array
  is not ahead by location in any band (10-20 mm: 7/4, p = 0.40; p = 0.52 in v3, 0.04 in v2 and 0.14
  in the earlier run: not robust).

Localization design (`scripts/g4_localization.py`, bounded)
* 24 locations (2 per depth x orientation stratum), focal dipoles and 10-mm patches at 80 and
  320 nAm, one event per location and condition (96 events), each in 4 s of independent null data
  shared by the arrays; truth as in detection (3-layer BEM, exact geometry).
* Inverse model with bounded mismatch: 1-layer BEM (inner skull) and a coregistration error of
  2 mm and 2 deg, drawn 8 times (translation direction, rotation axis) and shared by all arrays;
  location i uses draw i mod 8 in every condition (v2; v1 cycled draws by event, so each condition
  saw only 2 of them), median displacement at the true source 2.9 mm; noise covariance from 5 min
  of independent null data (Ledoit-Wolf).
* MNE/dSPM: the study's own implementation of the minimum-norm estimator (`localization.MNEInverse`:
  depth weighting with exponent 0.8 on the fixed-orientation whitened gains, no MNE-style limit on the
  weights, regularisation at SNR 3, dSPM normalisation by each source's noise standard deviation)
  and, from v4, MNE's own on the same events and grid (`mne.minimum_norm.make_inverse_operator` with
  fixed orientation on the discrete source space, depth 0.8 with MNE's default limit and its
  gradiometer-based depth weighting, then `apply_inverse`, dSPM, lambda^2 = 1/9), so the deviation is
  measured, not assumed. On a 5-mm Poisson-disk grid of 3,821 usable
  vertices that excludes the true source vertices; ECD with MNE's `fit_dipole` (same BEM and
  transform, at least 5 mm inside the inner skull) at the spike peak. Errors are measured on the
  MRI through the analyst's (perturbed) transform; the ECD error in the sensor frame is kept as a
  decomposition. Each event also passes through the practical detector (1 false event per
  minute, thresholds from 10 min of independent null data; from v4 their rate is checked on 10 min
  of independent held-out null data drawn from its own random stream, so the events are those the
  primary study would draw).
  The inverse is evaluated at the true peak sample (oracle timing), so the joint
  detection-and-localization success is conditional on the event time being known; localization at
  the detected event's time would add the detector's timing error.
* No goodness-of-fit cut: MNE computes GOF on whitened data, where noise adds about one unit per
  channel, so at equal SNR it is higher for arrays with fewer channels (median for detected 320-nAm
  focal events: 65 % matched OPM, 51 % dense OPM, 39 % Neuromag). GOF, chi2/dof and the 95 %
  confidence volume are descriptive. Failed fits: `mne.fit_dipole` returned a dipole for every event
  (96 events x 5 views in each of the nine anatomies; a failure would have stopped the run), so
  failures are counted under declared criteria from the stored per-event tables
  (`scripts/study_g4_fit_failures.py`, `results/g4/G4_fit_failures_report.md`; the three primary
  arrays): a dipole or dSPM peak more than 30 mm from the true source (gross), and a dipole
  confidence volume above 10 cm^3. Over the nine anatomies (864 events per array) gross dipole errors
  number 332 (Neuromag), 325 (matched OPM) and 310 (dense OPM), 27, 38 and 44 of them among detected
  events; gross dSPM errors 411, 327 and 329 (detected 77, 36 and 39); 87 % of all gross errors belong
  to undetected events and 73 % to 80-nAm sources. The medians reported below include these events
  ("all") or exclude the undetected ones ("detected"); 111 of the adult's 288 dipole errors exceed
  30 mm, 93 of them at undetected events and 81 of those at 80 nAm.
* Paired OPM-minus-Neuromag comparisons on identical events (one per location): median error
  difference with a bootstrap CI and Wilcoxon signed-rank p; exact McNemar p for joint detection +
  localization within 10 mm; 20 comparisons per OPM array against Neuromag combined (16 before the
  MNE dSPM was added), uncorrected. Secondary (v4): Neuromag also with its magnetometers or its
  gradiometers alone, each with its own noise covariance, detector, thresholds and inverse (the data
  are the combined system's), compared with the OPM arrays in the same way (`paired_secondary`).
* Coverage limits of the epilepsy examples: the 72 detection locations are frontal-heavy (24 frontal,
  2 occipital, 6 insular) and are not analysed by region; one patch extent (10-mm radius) and one
  morphology family (three stretches of one spike-wave complex).

Localization results (v4 arrays, run at e53bea8, `results/g4/g4_localization_summary.json`;
Neuromag combined, matched OPM, dense OPM)
* Detection in this subset: 320-nAm focal 0.92, 0.96, 1.00; 320-nAm patches 0.83, 0.79, 0.88;
  80-nAm patches 0.04-0.12 (too weak: their "localization" is that of noise). The 1-per-minute
  thresholds (10 min of null data) give 1.1-1.6 false events per minute on 10 min of held-out null
  data (exact 95 % intervals about 0.6-2.6): with about 10 calibration events per threshold the
  operating point is approximate, as for the detection study.
* ECD, detected events: 320-nAm focal 3.7, 5.8, 4.2 mm on the MRI (2.5, 4.2, 1.8 mm in the sensor
  frame; the coregistration error alone displaces the source by a median 2.9 mm); 320-nAm patches
  7.2, 8.1, 8.0 mm; 80-nAm focal 3.7, 6.2, 5.1 mm. chi2/dof 0.8-1.1: the head-model and
  coregistration mismatch is small against the noise at these SNRs.
* dSPM peak, detected events (the study's estimator | MNE's `mne.minimum_norm`): 320-nAm focal
  13.6 | 13.6, 10.2 | 10.0, 13.0 | 11.8 mm; 320-nAm patches 16.7 | 14.4, 12.9 | 12.5, 12.5 | 9.3 mm. The
  two implementations differ by a few millimetres in either direction, MNE's mostly lower (its depth
  weighting is limited and uses the gradiometers). Support recovery of the 320-nAm patches (share of
  the N strongest grid sources that are members of the patch, N being its number of grid sources;
  chance ~N/3,821): medians 0.0, 0.07, 0.07.
* Detected and localized within 10 mm (dSPM | ECD): 320-nAm focal 0.21 | 0.67, 0.46 | 0.58,
  0.33 | 0.67; 320-nAm patches 0.00 | 0.58, 0.29 | 0.50, 0.33 | 0.58.
* Paired (24 locations; 20 comparisons per OPM array against Neuromag combined, uncorrected p):
  for 320-nAm patches the dense array's dSPM error is 5.6 mm smaller (p = 0.054; with MNE's
  dSPM 1.0 mm, p = 0.22), the matched array's 2.4 mm (p = 0.067; MNE 2.3 mm, p = 0.073), and both
  are detected and dSPM-localized within 10 mm more often (8/0 and 7/0 locations, p = 0.008 and
  0.016), which does not survive a correction over the 20 comparisons (threshold 0.0025). The dipole
  errors do not differ (|difference| <= 1.8 mm, p >= 0.40). Earlier runs: v3 0.0 mm (p = 0.50 and
  0.81), v2 -4.8 and -5.3 mm (p = 0.14 and 0.046), an earlier v2 run -4.8 and -3.7 mm: a dSPM
  benefit of the OPM arrays for extended sources appears in some runs and not others, and its size
  depends on the inverse implementation. Against Neuromag's gradiometers alone (secondary) the OPM
  arrays' dSPM errors are smaller in several conditions (e.g. 320-nAm focal: matched -4.6 mm,
  p = 0.0003; dense -2.3 mm, p = 0.013); against the magnetometers alone 2 of the 40 comparisons have p < 0.05.
* With 24 locations the study bounds rather than resolves localization differences. In this model
  the dipole fit is limited by the coregistration error (about 3 mm) for every array; a dSPM
  benefit of on-scalp sensors for extended sources is not established.

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
  / '18mo' / '12mo')`; averages of many MRIs of each age, built from the Neurodevelopmental MRI
  Database of Richards et al. 2016) are used in their native dimensions with their own 3-layer
  BEMs, dense head surfaces, oct-6 source spaces (whose full white surfaces are the
  full-resolution cortices), aparc labels and fiducials (head frames from them). The templates are
  distributed publicly by their authors (O'Reilly et al. 2021; J. E. Richards, who created the
  database, is a co-author) under LGPL-2.1 through MNE-Python's `fetch_infant_template`; only
  results derived from them are published, citing O'Reilly et al. (2021) and Richards et al.
  (2016) wherever they appear. The 2-year template is the primary pediatric
  anatomy (the 2-year size control is scaled to its head circumference); the 18- and 12-month
  templates (added 2026-10-01) show how the results move between averages of one database, not
  between children. They are templates, not individual children: results on them are conditional
  simulations, not population estimates, and anatomical variability is not assessed. The 12-month
  template's skull is 0.25 mm thick at its thinnest (as distributed; 23 of 2,562 inner-skull
  vertices within 1 mm of the outer skull); the 1-layer BEM check does not depend on it.
* School-aged children (added 2026-10-03; D-G3-ANAT): three typically developing children of OpenNeuro
  ds005234 v2.2.0 (doi:10.18112/openneuro.ds005234.v2.2.0; Fadeev et al. 2024, J Neurodev Disord
  16(1):67), child A (sub-Z213, 7.8 years), B (sub-Z209, 8.3) and C (sub-Z226, 8.7),
  individual MRIs processed with FreeSurfer by the dataset's authors. Used: their own white and
  sphere surfaces, aparc annotations and dense MRI scalp (lh.seghead), which share each child's
  volume information. The snapshot's file tree shifts each subject's watershed BEM and talairach.xfm
  into the preceding subject's folder (the S3 objects behind them carry the name of the subject
  they belong to; the BEMs
  were identified by the same volume information as the child's surfaces). The watershed BEMs put
  the inner skull just below the scalp (a median 0.8-2.5 mm over the upper head, against 9.7 mm in
  the adult and 5.6-8.5 mm in the templates; `scripts/study_school_anatomy.py`,
  `results/g3b/school_anatomy_checks.json`), with an outer skull crossing it. The skull is therefore
  modelled (A-BEM-CHILD): the inner skull moved inward to 8 mm below the scalp where it is
  shallower, no vertex closer than 2 mm to a white-surface vertex (the surface between vertices comes
  within 0.9-1.3 mm of the white surface; no white vertex lies outside it), the outer skull halfway to
  the scalp; over the upper head the modelled inner skull lies a median 7.5, 5.7 and 7.8 mm below the
  scalp (child B's cortex lies close to its scalp). The head surface is the watershed outer skin put on the scalp (A-BEM-CONFORM). On the
  adult, whose segmented skull is known, the same procedure applied to a degraded inner skull changes
  the G2 headline by at most 0.003 (within 0.01 in every depth band) for depths of 6-10 mm, although
  at 8 mm its inner skull lies a median 1.7 mm from the segmented one (`scripts/study_child_bem.py`).
  No fiducials come with the data. Fiducials from MNI coordinates (fsaverage's, mapped through the
  talairach.xfm) are not used: on the adult, whose fiducials are digitised, they land 10-27 mm from
  them (head frame 10.6 deg off), and only child A's own talairach.xfm was obtained (with it child A's
  cortex centroid lies 3.9 mm from the adult's in MNI space; the files stored in the children's own
  folders, which belong to other subjects, put the three children's centroids 17-27 mm away). Instead
  the adult's digitised fiducials are transferred by a similarity fit of the cortices and moved to the
  nearest scalp vertex (A-G3-FID; on the templates this reproduces their own fiducials to 1.6-14.4 mm
  and their head frame to 3.2-5.9 deg, where the same fit between the scalps misplaces the nasion by
  34-45 mm). Prepared by `scripts/prepare_school_subjects.py` in the templates' file layout
  (`results/g3b/school_subjects_preparation.json`); FreeSurfer files are read without nibabel
  (`opmsquid.fsio`, checked against MNE's files of the sample subject). They are individuals, three
  of one dataset, not a population; like the templates they are compared with the adult by parcel and
  stratum, with their own bootstrap stream. No cortex maps are drawn for them (their inflated surfaces
  were not obtained).
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
  device x and y, pitch it by +-10 deg, roll it by +-5 deg or turn it by +-10 deg (yaw, v4) about the
  head origin, then raise it
  to the same top contact where there is room (an adult head shifted towards the helmet wall
  stays where it is). After a later check two variants were added: 'x-centred' shifts the
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
  projected sites are at least 21.4 mm apart). The 204-site channel-budget control has no
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
  (for the scaled controls vertex-wise, for the templates the adult-area-weighted median over parcels
  of the parcel-level difference, which differs from the difference of the two pooled medians)
  reported together. D does not depend on the moment. Peak-channel and mean-power SNR (dB) are
  secondary metrics. Homologous comparison: vertex-wise for the scaled controls; for the templates
  (no vertex correspondence) per Desikan-Killiany parcel and per declared depth (native mm below
  the scalp) and orientation stratum, from area-weighted medians (target weights = usable cortical
  area of each target's nearest-target cell). Strata or parcels with fewer than 10 targets in
  either anatomy are reported as sparse. Intervals: bootstrap over parcels (one anatomy, or each
  anatomy independently for between-anatomy strata); they contain no between-subject variability.
* Extended sources (A-G3-PATCH): fixed-total 10-nAm geodesic patches of 5, 10 and 20 mm (native)
  around 300 targets (every third of them for 20 mm; the same vertices in the adult and the scaled
  controls, random targets on each template), for the primary placement and the dense OPM array.
  Absolute detectability of the patches also at a fixed current density of 0.5 nAm/mm^2 (moment =
  density x patch area; human neocortex 0.16-0.77 nAm/mm^2, Murakami & Okada 2015).
* Usefulness (A-G3-USEFUL; operational, not a clinical standard): a source is usable by a system
  when d reaches 5 for a reference moment (20, 50, 100, 200 nAm; maps at 100 nAm); shares of the
  usable cortical area where OPM, SQUID, both or neither are usable, and the moment needed for
  d = 5 by depth.
* The medial wall (FreeSurfer 'unknown': 558 adult targets; 575, 631 and 496 on the 24-, 18- and
  12-month templates) is not cortex and is left out of every G3B summary (it made up most of the
  60-90 mm strata).

G3B results (`results/g3b/g3b_summary.json`, `G3B_report.md`; computed at 71de176 with all nine anatomies and
the v4 arrays: equal 7-mm standoff, A-BEM-CONFORM; the six anatomies of the v4 pass, computed at e53bea8 and
redrawn at 9512f1b, reproduce exactly, 97,048 values; the school-aged children have their own bootstrap
stream; summaries, report and figures redrawn at 6ebac95 to mark an infeasible placement and the notes
at ff6e25c and d3bb1a9, every number unchanged). For the six, against v3 (19a8fd2) the adult's D rose from +0.85 to +1.00 dB and the scaled
controls' by 0.04-0.05 dB, while the templates, whose head surfaces already lay on their scalps, did not
change; every Delta fell by 0.11-0.27 dB and the conclusions did not change, below. Dense OPM vs
Neuromag, intrinsic + brain noise, primary placement unless stated; dB of detectability;
area-weighted medians without the medial wall, parcel-bootstrap 95 % intervals. "Template" alone
means the 2-year template; the 18- and 12-month templates are named.
* Heads: occipitofrontal circumference 586 mm (adult), 525 (school-age size, x0.895), 495 (2-year
  size, x0.844) and 495 mm (template); usable cortex 1,878, 1,470, 1,290 and 1,062 cm^2. Refitted
  OPM arrays: dense 208, 174, 155 and 151 sites (23.1-25.9 per 100 cm^2 of covered scalp, which
  shrinks from 901 to 584 cm^2), matched 98, 89, 90 and 83 (built from the Neuromag sites at the
  primary placement, which for the adult is 5.5 mm above G2's measured position; as many sites pass
  the site rules there as in G2, 98). Top contact raises the heads by 5.5,
  22, 22 and 28 mm; the median magnetometer-to-scalp gap is then 28.4, 34.1, 39.1 and 38.5 mm
  (centred 29.8, 41.3, 46.9, 47.3 mm). Every placement of these four (17 per anatomy with the counterfactual
  helmets) is feasible;
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
* Link to G2: at the adult's measured position the dense/combined ratio is 1.144x as G2 reports it
  (unweighted, all targets), 1.148x without the medial wall and +1.27 dB area-weighted; at top
  contact the adult's D is +1.00 dB [+0.82, +1.17] (1.12x).
* D_child, D_adult, Delta (vs Neuromag combined): school-age size +1.50 vs +1.00 dB, Delta +0.44
  [+0.34, +0.56] dB; 2-year size +2.13 vs +1.00, Delta +1.07 [+0.91, +1.18]; 2-year template +1.84
  vs +1.00, Delta +0.73 [+0.46, +0.98] (66 parcels, Delta > 0 in 93 % of their area; the frontal
  poles are sparse). Against the gradiometers alone Delta is +0.65, +1.49 and +1.13 dB, against the
  magnetometers alone +0.41, +0.96 and +0.82 dB. With the external projection +0.44, +1.09 and
  +0.42 dB; for the matched-site OPM array +0.34, +1.10 and +0.68 dB (the gain is not a coverage
  effect of the dense array); for 5-, 10- and 20-mm patches +0.58/+0.54/+0.44, +1.30/+1.22/+0.99 and
  +0.56/+0.58/+0.72 dB.
  Secondary metrics: peak-channel SNR (D_adult -1.11 dB, Neuromag's best channel ahead) Delta
  +0.69, +1.66, +1.52 dB; mean-power SNR (D_adult +0.34 dB) Delta +0.39, +0.74, +0.83 dB.
  18- and 12-month templates: D_child +1.89 and +2.04 dB, Delta +0.88 [+0.54, +1.05] and +0.96
  [+0.61, +1.51] dB (66 and 65 parcels, Delta > 0 in 96 and 97 % of their area); against the
  gradiometers +1.26 and +1.33, the magnetometers +0.72 and +0.98 dB; projected +0.57 and +0.83 dB;
  matched-site array +0.61 and +0.76 dB; 5-, 10- and 20-mm patches +1.05/+0.81/+0.89 and
  +1.05/+1.02/+0.91 dB;
  peak-channel SNR +1.32 and +1.60, mean-power SNR +0.92 and +1.08 dB. Across the three templates
  Delta is +0.73, +0.88 and +0.96 dB (24, 18, 12 months; overlapping intervals), but reweighted to
  the adult's depth mix the pooled differences are +0.32, +0.53 and +0.29 dB (below): the
  12-month template's top rank reflects its shallow cortex (the raw ordering also follows head
  circumference, 495, 491 and 469 mm, so depth mix and size are not separated here). The templates'
  shallower cortex (42 % of the 2-year template's area at 10-20 mm depth vs 21 % of the adult's) may
  partly be an artefact of template averaging, which smooths sulci; the depth-reweighted differences
  are the more conservative statement.
* School-aged children (individual MRIs; A, B, C at 7.8, 8.3 and 8.7 years): occipitofrontal
  circumference 520, 486 and 535 mm (school-age size control 525 mm, 2-year size control 495 mm), but
  usable cortex 1,852, 1,666 and 1,869 cm^2, close to the adult's 1,878 cm^2 (school-age size control
  1,470 cm^2);
  dense 155, 153 and 167 sites, matched 90, 85 and 91; top contact raises them by 17, 26.5 and 21 mm
  (median gap 35.2, 38.4 and 32.6 mm; centred 40.2, 47.3 and 41.3 mm). D_child +1.34, +1.51 and +1.19
  dB, Delta +0.30 [+0.11, +0.38], +0.41 [+0.28, +0.54] and +0.17 [-0.07, +0.34] dB (66 parcels each;
  Delta > 0 in 71, 84 and 59 % of their area): at a similar head circumference children A and C gain
  less than the school-age size control (+0.44 [+0.34, +0.56] dB; the intervals overlap or touch),
  child B, with a smaller head, about as much. Against the gradiometers +0.18, +0.34 and +0.11, the magnetometers +0.33, +0.50 and
  +0.25 dB; with the external projection +0.03, +0.26 and -0.19 dB; matched-site array +0.27, +0.08 and
  +0.17 dB; 5-, 10- and 20-mm patches +0.20/+0.31/+0.50, +0.49/+0.57/+0.35 and +0.06/+0.06/+0.35 dB;
  peak-channel SNR +0.64, +0.98 and +0.24, mean-power SNR +0.16, +0.23 and +0.17 dB. Their cortex is
  shallower than the adult's (42, 47 and 35 % of their area at 10-20 mm vs 21 %; median depth 21.6, 20.0
  and 23.4 vs 26.2 mm): reweighted to the adult's depth mix the pooled difference falls from +0.34, +0.52
  and +0.19 to -0.04, +0.05 and -0.05 dB, and within depth strata Delta is near zero or negative down to
  30 mm (-0.36 to +0.08 dB; intervals excluding 0 for child A at 20-30 mm, child B at 15-25 mm and child
  C at 20-25 mm) and positive deeper (30-60 mm: +0.01 to +1.77 dB). Mechanism: with the background
  fixed per unit area, their nearly adult-sized cortex keeps Neuromag's brain noise near the adult's
  (magnetometers 194, 168 and 174 fT vs 202 fT; gradiometers 40, 34 and 36 vs 41 fT/cm) and raises the
  closer dense OPM's (711, 751 and 617 fT vs 494 fT) along with its peak signal (293, 326 and 240 fT vs
  208 fT); absolute detectability of a 10-nAm dipole: dense OPM -1.55, -1.41 and -2.11 dB (adult -1.60),
  Neuromag combined -2.91, -2.70 and -3.26 dB (adult -2.67). What drives Delta: centred +0.76, +1.57 and
  +1.00 dB; laterally centred, then top contact +0.29, +0.36 and +0.09; true 18-mm contact +0.31, +0.40
  and +0.17; back +0.44, +0.94 and +0.87 dB; in the counterfactual helmet -0.74 [-0.86, -0.60], -0.23
  [-0.49, -0.05] and +0.56 [+0.37, +0.78] dB (factors 0.897, 0.884 and 0.973, raised from the
  head-circumference ratios 0.887, 0.829 and 0.913 until no magnetometer is within 18 mm of the scalp:
  child C's helmet barely shrinks), about the laterally centred head -0.76 [-0.84, -0.68], -0.59
  [-0.74, -0.42] and +0.04 [-0.20, +0.20] dB (factors 0.897, 0.864 and 0.943). At an equal channel count Delta would be +0.66, +0.78 and +0.43 dB. Top contact
  ranks 4th, 6th and 3rd from the lowest of the 12 source-blind placements (child C: of the 11 feasible;
  its x-5mm placement brings a magnetometer to 16.9 mm from the scalp, inside the 18-mm Dewar spacing,
  and is marked infeasible) (D_child 1.25-1.50, 1.40-1.96 and 1.05-1.91
  dB over the others; against the family medians the differences of the medians are 0.12, 0.17 and 0.05
  dB smaller). Sensitivity (difference of the medians): OPM noise 7-30 fT/sqrt(Hz) -0.15 to +0.60, -0.03
  to +0.77 and -0.10 to +0.32 dB (negative at 7 fT/sqrt(Hz)), background x0.5 / x2 +0.46 / +0.20, +0.67
  / +0.34 and +0.25 / +0.11 dB, 1-layer BEM +0.35, +0.48 and +0.19 dB. At 0.5 nAm/mm^2, 10-mm patches
  reach d >= 5 at 40 / 51, 43 / 53 and 36 / 50 % of the centres (Neuromag / dense OPM; adult 43 / 54 %).
* Equal standoff (v4): every anatomy's head surface lies on its MRI scalp (A-BEM-CONFORM), so the
  OPM arrays of all nine heads sit at a median 6.99-7.01 mm above their MRI scalp (sites moved out by
  the clearance rule: adult 26 of 208, scaled controls 18 and 17, templates 3-5, school-aged children
  6, 3 and 4). In v3 the adult's
  outer skin lay about 1 mm outside its scalp, so its sensors sat about 0.8 mm farther out than the
  children's and its axes were tilted (sections 2.2, 3). The v3 sensitivity that moved the children's
  arrays outward by the adult's excess predicted Delta +0.47, +1.12, +0.77, +0.81 and +1.07 dB and a
  counterfactual Delta about the laterally centred head of -0.34 to -0.06 dB; v4 gives +0.44, +1.07, +0.73, +0.88 and +0.96 dB and
  -0.38 to -0.16 dB (school-age size, 2-year size, 24, 18, 12 months): Delta within 0.11 dB of the
  prediction, and a counterfactual Delta about the laterally centred head that is now negative in each
  of these five heads (the templates' intervals include 0). The adult's D changed through its
  standoff, its sensitive axes and its lower occipital sites together (section 8).
* By depth (template vs adult, combined, native depth strata): Delta +0.17 to +0.67 dB down to 30 mm
  (the 10-20 mm intervals exclude 0, the 20-30 mm ones include it), +0.29 [+0.12, +0.55] at 30-40 mm, +0.24 at 40-50 mm, +0.85 at 50-60 mm and
  +1.62 at 60-90 mm (83 template vs 26 adult targets, 79 of them isthmus cingulate). Within strata the
  template's Delta is largest deep and for radial sources, but those strata hold little of its area
  (2.9 % of it is deeper than 50 mm); most of its pooled gain reflects its shallower cortex (42 % of its area at
  10-20 mm vs 21 % in the adult, area-weighted median depth 21.8 vs 26.2 mm): reweighted to the
  adult's depth mix, the target-level difference of the medians falls from +0.85 to +0.32 dB. For the
  scaled controls the homologous (vertex-wise) Delta is high near the surface and falls with depth
  to 40-60 mm: 2-year size +2.40 [+1.78, +2.54] dB at an adult depth of 10-15 mm, +1.57 at 20-25 mm,
  +0.64 at 30-40 mm, +0.35 at 40-50 mm; school-age size +1.05, +0.76, +0.26 and +0.13 dB; it rises
  a little in the deepest stratum, which holds 26 targets (60-90 mm: +0.71 and +0.33 dB, intervals
  include 0). By
  orientation: the scaled controls gain somewhat more for tangential than for radial sources (+0.41 /
  +0.58 dB; +0.93 / +1.29 dB; overlapping intervals), the template far more for radial sources (0-30 deg +1.94 [+1.33,
  +2.50], 30-60 deg +0.93, 60-90 deg +0.59 dB; radial sources are 16 % of its targets by count). These are
  pooled over depth; at matched depth radial sources still gain more down to 40 mm (+0.76, +0.65,
  +0.63 vs tangential +0.22, +0.30, +0.11 dB at 0-15, 15-25, 25-40 mm) but not deeper (+0.08 vs
  +0.61 dB at 40-90 mm). The 18- and 12-month templates repeat this, except that the 12-month template's Delta is largest at
  10-15 mm (`template_depth_checks`):
  within strata Delta is +0.22 to +1.15 dB at 10-50 mm (every interval but the 12-month
  template's at 20-25 mm excludes 0; 12 months at 10-15 mm +1.15 [+0.45, +1.77]) and +0.78 to +1.52 dB deeper; reweighted to the adult's
  depth mix their pooled difference falls from +0.90 to +0.53 and from +1.05 to +0.29 dB (36 and 46 %
  of their area at 10-20 mm; median depth 23.2 and 20.6 mm); at matched depth radial sources gain
  more than tangential ones down to 40 mm (18 months +1.19, +1.32, +0.88 vs +0.38, +0.49, +0.20 dB;
  12 months +1.49, +0.70, +0.64 vs +0.52, +0.41, +0.10 dB) but not deeper (+0.11 and
  +0.08 vs +0.52 and +0.57 dB at 40-90 mm).
  Regionally the template's Delta is asymmetric (left inferior temporal, pars orbitalis,
  entorhinal, middle temporal and fusiform +2.0 to +3.0 dB; right orbitofrontal -0.9 to -1.0 dB): at top contact under the
  adult's measured pose its left frontal, temporal and parietal gaps are about 8 mm wider than the
  right ones (occipital 1.7 mm); the x-5mm variant, which roughly re-centres it, changes the median
  D by -0.1 dB.
* What drives Delta (each child placement against the adult at the same rule): with the child left
  centred (ear line where the adult's was) Delta is +1.14, +1.68 and +1.83 dB; at top contact +0.44,
  +1.07 and +0.73 dB; laterally centred, then top contact, +0.44, +1.06 and +0.71 dB; at true 18-mm
  contact +0.44, +1.09 and +0.77 dB; back contact +0.85, +1.21 and +1.39 dB. In the counterfactual
  helmet scaled with the head it is -0.22 [-0.25, -0.18], -0.37 [-0.43, -0.33] and +0.24 [+0.13,
  +0.56] dB, and about the laterally centred head -0.21, -0.38 and -0.30 [-0.42, +0.02] dB: the
  template's positive residual came from its lateral offset. For the 18- and 12-month templates:
  centred +2.08 and +2.50, top +0.88 and +0.96, laterally centred then top +0.75 and +0.91, 18-mm
  contact +0.85 and +0.85, back +1.49 and +1.68 dB; counterfactual +0.18 [+0.11, +0.47] and -0.02
  [-0.16, +0.11] dB, about the laterally centred head -0.16 [-0.39, +0.03] and -0.18 [-0.37, +0.06]
  dB. The relative gain of the head-adaptive array in the templates and size controls therefore comes
  from the fixed helmet's fit; with a helmet that fits as the adult's does, it reverses (about the
  laterally centred head -0.21 and -0.38 dB for the scaled controls; -0.30, -0.16 and -0.18 dB for the
  templates, whose intervals include 0; school-aged children A and B -0.76 and -0.59 dB, child C level
  at +0.04 [-0.20, +0.20] dB).
* Mechanism (templates and size controls; for the school-aged children, whose cortex is nearly
  adult-sized, see their bullet above). Both systems' detectability rises in these smaller heads, by
  different routes:
  vertex-wise from the adult, the dense OPM gains +1.10 and +1.67 dB in the two scaled controls,
  Neuromag combined +0.65 and +0.54 dB (gradiometers +0.45 and +0.10, magnetometers +0.69 and
  +0.64 dB). The on-scalp OPM sees more signal from a cortex closer in absolute terms (median peak
  field of a 10-nAm dipole 208 fT in the adult, 241, 260 and 271 fT in the children; 255 and 293
  fT at 18 and 12 months) while its brain noise stays near 450-540 fT. In the fixed helmet the
  cortex is farther from the SQUIDs (median target-to-magnetometer distance at 15-20 mm depth 48
  mm adult, 50-57 mm children; OPM 26-27 mm in every head) and, with the background fixed per unit area, the smaller cortex carries
  less background power, so their brain noise falls (magnetometers 202 -> 150, 125 and 120 fT, 110
  and 98 fT at 18 and 12 months; gradiometers 41 -> 30, 23, 23, 21 and 18 fT/cm against 21 fT/cm
  intrinsic) more than their signal (magnetometer peak 70 -> 66, 54, 65, 63 and 67 fT). Absolute
  detectability of a 10-nAm dipole (median dB): dense OPM -1.60 (adult), -0.61, -0.20, +0.92
  (template), +1.07 (18 months) and +1.86 (12 months); Neuromag combined -2.67, -2.15, -2.27,
  -1.19, -1.03 and -0.62. In the counterfactual helmet about the laterally centred head the SQUID
  gains at least as much as the OPM (also in children A and B; child C +0.04 [-0.20, +0.20] dB): the
  counterfactual Delta depends on the comparator (gradiometers -0.37
  and -0.60 dB, magnetometers -0.16 and -0.28 dB for the scaled controls), which points to a
  SQUID-side component (signal against fixed intrinsic noise, largest for the gradiometers, whose
  brain noise is only about twice their intrinsic noise) rather than to the OPM's fixed standoff
  alone; the noise composition in the counterfactual helmet is not reported.
* Channel count: the children's dense arrays have fewer sites (174, 155, 151). The adult's dense
  array subsampled (farthest-point) to those counts has D +0.73, +0.59 and +0.54 dB (full array
  +1.00), so at an
  equal channel count Delta would be +0.63, +1.41 and +1.23 dB (18- and 12-month templates, 157 and
  144 sites: adult subsampled +0.61 and +0.46, Delta +1.37 and +1.43 dB): the smaller site count of
  the head-adaptive array works against it, and the headline Delta includes that loss.
* Extended sources: the 5- and 10-mm patches' D is within 0.15 dB of the focal D at the same centres
  in every head; for the 20-mm patches it is 0.1-0.7 dB lower in the adult, the size controls, the
  templates and child A (adult +0.83 vs +1.05 dB), 0.08 dB lower in child B and 0.23 dB higher in
  child C, and Delta is +0.44, +0.99, +0.72, +0.89 and +0.91 dB (children A-C +0.50, +0.35 and
  +0.35 dB). At a fixed current density of 0.5 nAm/mm^2 (median moments about 35, 140 and 570 nAm)
  the share of patch centres with d >= 5 is, Neuromag / dense OPM over the nine heads: 5-mm patches
  0-8 % / 6-32 %, 10-mm 31-77 % / 50-85 % (adult 43 / 54 %), 20-mm 85-100 % / 91-100 %; both are
  higher on the templates, whose cortex is shallower.
* Placements: the source-blind variants (+-5 mm, pitch +-10 deg, roll +-5 deg, yaw +-10 deg, back
  contact) give
  D_child 1.46-2.56 (school-age size), 1.69-3.04 (2-year size), 1.75-2.51 (template), 1.80-2.48
  (18 months) and 1.98-2.79 dB (12 months), adult 0.98-1.43 dB; top contact (the primary) ranks in
  the middle of the family for the templates and size controls (5th-7th of 12; children A-C 4th, 6th
  and 3rd, above) but second-lowest for the adult, so against each anatomy's family median the
  child-minus-adult difference of the medians is 0.14-0.20 dB smaller (+0.33, +1.00, +0.67, +0.70 and
  +0.88 dB vs +0.51, +1.14, +0.85, +0.90 and +1.05 dB at top contact; children A-C 0.05-0.17 dB).
  Among the placements, true 18-mm contact gives the lowest D for the adult, the school-age size
  control, the three templates and children A and B (0.87, 1.38, 1.70, 1.75, 1.85, 1.24 and 1.40 dB),
  but not for the 2-year size control (2.03 dB; pitch +10 deg gives 1.69) or child C (1.08 dB; pitch
  +10 deg gives 1.05). Regions:
  raising the head (centred -> top) lowers D in every lobe, most in the parietal and frontal
  lobes (template frontal +3.71 -> +2.02, parietal +3.02 -> +1.39 dB); back contact favours the
  SQUID at the occiput (template occipital +1.25 dB vs +2.08 at top) and disfavours it frontally
  (+3.89 vs +2.02 dB).
* Sensitivity (difference of the median D, child minus adult, with the same variant applied to both;
  a check on the medians, not the paired Delta estimator; vs combined): OPM noise 7-30 fT/sqrt(Hz)
  +0.42 to +0.52 (school-age size), +1.03 to +1.14 (2-year size), +0.66 to +0.85 dB (template);
  background variance x0.5 / x2 +0.53 / +0.48, +1.18 / +1.07, +0.89 / +0.78 dB (the variants scale
  the adult too; D_child itself moves by at most 0.06 dB, template 1.86 / 1.79 vs 1.84 dB; in the
  school-aged children by up to 0.16 dB, child B 1.64 / 1.35 vs 1.51 dB); 1-layer
  BEM +0.51, +1.13, +0.83 dB. 18- and 12-month templates: OPM noise +0.66 to +0.94 and +0.82 to
  +1.05 dB, background x0.5 / x2 +0.89 / +0.89 and +1.08 / +0.97 dB, 1-layer BEM +0.84 and +1.06 dB.
  At 30 fT/sqrt(Hz) the adult's D is -0.05 dB (Neuromag slightly ahead) and the children's +0.37, +0.97,
  +0.61, +0.61 and +0.77 dB (children A-C +0.55, +0.72 and +0.27 dB).
* Usefulness (d >= 5) at 100 nAm, share of usable cortical area, both / OPM only / SQUID only /
  neither: adult 0.66 / 0.02 / 0.00 / 0.32, school-age size 0.69 / 0.04 / 0.00 / 0.28, 2-year size
  0.69 / 0.05 / 0.00 / 0.26, template 0.75 / 0.05 / 0.00 / 0.20, 18 months 0.76 / 0.05 / 0.00 / 0.19,
  12 months 0.79 / 0.05 / 0.00 / 0.16, children A-C 0.65 / 0.04 / 0.01 / 0.30, 0.65 / 0.04 / 0.01 / 0.29
  and 0.63 / 0.04 / 0.01 / 0.32; at 50 nAm OPM only 0.07, 0.10, 0.14, 0.12, 0.12, 0.13, 0.10, 0.10 and
  0.09. Over the reference moments (20-200 nAm) the SQUID alone is usable on at most 0.14 % of the area
  in the adult, the size controls and the templates (2-year template at 100 nAm) and on 0.75-1.18 % in
  the school-aged children (at 100 nAm); the cortex neither reaches is mostly deep and medial.
* Limits: one adult, three templates of one database, three school-aged children of one dataset (their
  skull modelled, their fiducials transferred) and two scaled copies; variability between three children
  is shown, not estimated; the templates are averages with smooth cortices; no age-specific background physiology (only the bounded x0.5/
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
* Anatomies: the 2-year (primary), 18- and 12-month templates (O'Reilly et al. 2021, built from
  the Neurodevelopmental MRI Database of Richards et al. 2016), the two size-only
  controls (scaled adults) and the three school-aged children with individual MRIs (A, B and C:
  7.8, 8.3 and 8.7 years; OpenNeuro ds005234 v2.2.0, Fadeev et al. 2024, with the modelled skull
  A-BEM-CHILD and the
  transferred fiducials A-G3-FID, section 10). Templates are averages and three children are not a
  population; the comparison with the adult mixes head size, anatomy and the fixed-helmet fit,
  which G3B separates for detectability.

Pediatric G4 results (`results/g4/G4_pediatric_report.md`, `g4_pediatric_comparison.json`; detection in
the adult and the five earlier children simulated at ed852b2 with the v4 arrays, their localization
at e53bea8, both studies of children A-C at 71de176, the comparison at 011cec1 (adding the children
left every earlier anatomy's entry unchanged); 18 locations per depth band in every anatomy; p-values
uncorrected, over 24 paired detection comparisons (3 comparators x 2 detectors x 4 bands) and 20
localization comparisons per OPM array and anatomy). The children's Neuromag is at the primary G3B
placement (top contact), the adult's at its measured position (the G4 adult study as frozen); top
contact would raise the adult's head by only 5.5 mm, against the children's 17-33.5 mm. Order below:
school-age size, 2-year size, 2-year template, 18 months, 12 months, then children A, B and C. The
templates' arrays did not change in v4, so their outcomes are essentially those of v3 (every
location count, p-value and S50 identical; the background scale, calibrated on the adult, changed
slightly, which moves some held-out rates by up to 0.05 per minute); the scaled controls' arrays did
(174 and 155 dense sites), and with their channel counts, so did their random streams.
* Held-out null: 0.4-1.55 false events per minute at the 1-per-minute thresholds (adult 0.65-1.3):
  thresholds calibrated on 20 min and checked on 20 min differ by up to this much, so operating
  points are approximate, and they are lopsided between arrays in either direction (dense OPM vs
  Neuromag combined: 1.05 vs 1.4, 1.2 vs 0.4, 1.3 vs 0.7, 1.1 vs 0.9, 0.85 vs 1.2, 0.75 vs 0.9, 0.7
  vs 0.8 and 0.75 vs 1.2 per minute). At a matched held-out rate of 1 per minute the superficial
  advantage remains (sensitivity for 40-nAm spikes at 10-30 mm, dense OPM vs Neuromag combined: 0.40
  vs 0.20, 0.31 vs 0.14, 0.43 vs 0.26, 0.50 vs 0.37, 0.44 vs 0.28, 0.23 vs 0.12, 0.17 vs 0.12 and
  0.31 vs 0.16; adult 0.36 vs 0.20). With every detector's threshold set on the held-out null to 1
  false event per minute (`scripts/study_g4_matched_rate.py`, `results/g4/G4_matched_rate_report.md`;
  in-sample for the held-out data) the paired 10-20 mm result is unchanged in every anatomy (the
  five earlier children: 13-16 locations favour the dense OPM, strength ratio 1.63, 1.56, 1.40, 1.33
  and 1.51, p <= 0.003; children A-C: 12/4, 13/3 and 14/2, ratio 1.29, 1.21 and 1.47, p = 0.015,
  0.0056 and 0.0031). Of the deeper differences below, the 12-month template's 30-45 mm (9/1, p =
  0.018) remains, and the adult's 20-30 and 45-70 mm (section 9); the 2-year size control's 30-45
  and 45-70 mm (8/3, p = 0.094; 7/2, p = 0.12) and the 18-month template's 45-70 mm (7/1, p = 0.15)
  do not; child C's 20-30 mm (10/2, p = 0.063, below) becomes 10/1 (p = 0.032). For the matched-site
  array, 10-20 mm: 10/3, 12/0, 9/7, 11/1, 11/4, 8/6, 10/5 and 8/4 (p = 0.064, 0.0005, 0.66, 0.004,
  0.03, 1.0, 0.59 and 0.42).
* Strength for 50 % detection, practical detector, Neuromag combined vs dense OPM, 10-20 mm: adult
  47 vs 34, then 53 vs 33, 61 vs 37, 43 vs 30, 36 vs 27, 43 vs 29, 61 vs 48, 65 vs 55 and 52 vs 36
  nAm; 45-70 mm adult 302 vs 254, then 299 vs 275, 273 vs 247, 240 vs 229, not reached vs 282 and
  254 vs 248 nAm; in children A-C Neuromag's is not reached and the dense OPM's only in child C (320
  nAm). The templates' point estimates are lower than the adult's at 10-30 mm and children A-C's
  higher at 10-45 mm, but the intervals overlap (e.g. Neuromag at 20-30 mm 71 [61-95], 62 [53-76]
  and 71 [57-98] nAm in the templates, 111 [78-160], 132 [107-181] and 95 [73-112] nAm in children
  A-C, vs 86 [66-116] nAm), and the sampled locations' median depths differ between anatomies by a
  few mm (10-20 mm: 13.6-17.8 mm, adult 16.9; 45-70 mm: 48.0-54.3 mm, adult 51.7): no difference
  between the anatomies is claimed.
* Paired, dense OPM vs Neuromag combined (locations favouring OPM / Neuromag; paired strength ratio
  Neuromag / OPM): at 10-20 mm every anatomy favours the OPM: adult 11/1 (1.36 [1.10-1.58]), then
  13/0 (1.60 [1.31-1.87]), 17/0 (1.65 [1.41-1.87]), 14/2 (1.44 [1.20-1.65]), 16/0 (1.33
  [1.16-1.62]) and 14/1 (1.46 [1.21-1.84]), p <= 0.004 in each; children A-C 11/4 (1.27
  [0.93-1.63], p = 0.022), 12/3 (1.18 [1.02-1.51], p = 0.0085) and 14/2 (1.43 [1.06-1.62], p =
  0.0037), none of the three surviving the correction over 24 comparisons (threshold 0.0021) that
  the five earlier children's results survive. Deeper, a location-level difference appears for the
  adult (20-30 mm 9/1, p = 0.014; 45-70 mm 9/0, p = 0.004; section 9), the 2-year size control
  (30-45 mm 8/2, p = 0.049, 1.16 [1.01-1.35]; 45-70 mm 8/1, p = 0.031, 1.11, interval open at both
  ends), the 12-month template at 30-45 mm (8/1, p = 0.031, 1.10 [1.03-1.23]) and the 18-month
  template at 45-70 mm (7/1, p = 0.047; Neuromag's S50 is not reached: ratio > 1.13, interval open
  at both ends); none survives a correction over 24 comparisons (threshold 0.0021; smallest p 0.004
  in the adult, 0.031 in the children), and elsewhere p >= 0.148 in the five earlier children
  (ratios 1.01-1.12) and p >= 0.0625 in children A-C (ratios 0.97-1.15; the smallest p child A's
  45-70 mm, 5/0, p = 0.0625, and child C's 20-30 mm, 10/2, p = 0.063). With the oracle every anatomy favours the
  OPM at 10-20 mm (11 to 17 locations favouring the OPM, at most 3 Neuromag), five of the nine at
  20-30 mm (none of children A-C: 9/3, 4/5 and 6/2) and four in a deeper band (the 2-year size
  control at 30-45 mm, 11/0, p = 0.001; the 2-year template at 45-70 mm, 11/0, p = 0.001; the
  18-month template at 30-45 and 45-70 mm, 10/3 and 12/2, p = 0.04 and 0.004; child B at 45-70 mm,
  8/0, p = 0.008; the 12-month template 9/1, 8/2 and 4/5 deeper than 20 mm). The matched-site array
  favours the OPM at 10-20 mm with the practical detector in the 2-year size control (14/0, p =
  0.0001) and the 18-month template (12/0, p = 0.0005), not established in the 12-month template
  (11/4, p = 0.051), the school-age control (8/3, p = 0.25), the 2-year template (10/5, p = 0.15) or
  children A-C (8/6, 11/5 and 8/4; p = 1.0, 0.27 and 0.54); unlike v3 (adult 1/9 and school-age
  control 0/7, p = 0.016), no anatomy has a matched-site deficit in its deepest band (adult 3/6,
  school-age control 4/5, children A-C 2/0, 3/2 and 4/3).
* So the detectability gains of G3B (Delta +0.44 to +1.07 dB in the templates and size controls,
  +0.17 to +0.41 dB in children A-C) are not resolved by the practical detector with 18 locations
  per band beyond the superficial band, whose advantage is present in the adult and every smaller
  head at a similar size (strength ratio 1.36 in the adult, 1.33-1.65 in the five earlier children
  and 1.18-1.43 in children A-C, where it is less certain; overlapping intervals); they appear in
  the oracle's 20-30 mm band (not in children A-C) and partly deeper.
* Localization (24 locations; Neuromag, matched, dense): ECD errors of detected events are similar
  in every anatomy (320-nAm focal: 3.7-7.6 mm; children A-C 4.4-6.7 mm). dSPM, all events, 320-nAm
  focal: the dense array is paired-closer than Neuromag in the school-age control (-7.0 mm, p =
  0.0002, surviving the within-anatomy correction below), the 2-year size control (-3.6 mm, p =
  0.02) and the 2-year and 18-month templates (-4.7 and -5.3 mm; p = 0.035 and 0.0033), not in the
  12-month template (0.0 mm, p = 0.38) or the adult (0.0 mm); in children A-C by -4.7, -2.9 and
  -3.0 mm (p = 0.18, 0.071 and 0.25); 320-nAm patches: school-age control -1.8 mm (p = 0.0049),
  2-year template -1.2 mm (p = 0.007), 12 months -2.0 mm (p = 0.0061), adult -5.6 mm (p = 0.054),
  children A-C -6.6, -0.1 and -4.6 mm (p = 0.22, 0.52 and 0.053). With MNE's own dSPM the dense
  array's 320-nAm focal differences are -5.0, -4.6, -1.8, 0.0 and -0.3 mm (p = 0.064, 0.054, 0.069,
  0.030 and 0.055; at 18 months 11 of the 14 untied locations favour the OPM): the same direction,
  smaller in the templates; in children A-C -5.5, -5.5 and -4.9 mm (p = 0.020, 0.009 and 0.012). The
  matched array in child A: dSPM -7.9 mm for 320-nAm and -10.3 mm for 80-nAm focal events (p =
  0.0046 and 0.0025). ECD, 80-nAm sources, 18 months: -3.0 mm (focal, p = 0.003) and -4.3 mm
  (patches, p = 0.002). Detected and localized within 10 mm (dSPM), 320-nAm focal spikes, Neuromag
  vs dense: adult 0.21 vs 0.33; children 0.21 vs 0.50, 0.12 vs 0.21, 0.29 vs 0.54, 0.29 vs 0.58,
  0.42 vs 0.58, 0.21 vs 0.38, 0.08 vs 0.33 and 0.33 vs 0.33. Within an anatomy and OPM array (20
  localization comparisons: the study's and MNE's dSPM errors, ECD errors and joint
  detection-and-localization success with either estimator, for focal and patch sources at 80 and
  320 nAm; the comparison file holds the 12 error comparisons, the localization summaries all 20)
  eight survive a Bonferroni correction (p < 0.0025), all favouring an OPM array: the matched
  array's dSPM of 80-nAm patches at 12 months (-33 mm, p = 0.0002; MNE's dSPM -27 mm, p = 0.0021;
  these weak patches are mostly not detected, so this compares noise-dominated estimates), the dense
  array's dSPM of 320-nAm focal events in the school-age control (-7.0 mm, p = 0.0002), MNE's dSPM
  of 80-nAm focal events at 12 months (matched -18 mm, p = 0.0004; dense -12 mm, p = 0.0017) and of
  80-nAm patches at 18 months (matched -17 mm, p = 0.0017), the matched dSPM of 320-nAm focal events
  at 18 months (-4.8 mm, p = 0.0019) and the dense ECD of 80-nAm patches at 18 months (-4.3 mm, p =
  0.0020). None survives in the 2-year template, the 2-year size control, children A-C (the
  smallest p there, 0.00252 for child A's matched dSPM of 80-nAm focal events, misses the threshold)
  or the adult; across the eight children and both arrays (320 comparisons, p < 0.00016) none (over
  the five earlier children's 200, p < 0.00025, the first two survived), and which comparisons
  survive varies between runs (v3 had four survivors under 16 comparisons, two of them among these;
  v2 seven). The direction is the more robust observation: over all nine anatomies 65 of the 360
  localization comparisons have p < 0.05 (uncorrected; about 18 would be expected by chance if they
  were independent, which they are not), 64 of them in favour of an OPM array (children A-C: 13 of
  their 120, all favouring an OPM array); 28 of the 65 concern weak 80-nAm sources, mostly undetected, 21
  are the study's dSPM errors over all events and 27 MNE's (without MNE's dSPM, as in v3: 38 of 288,
  37 favouring an OPM array).
* The pediatric epilepsy examples use the same framework as the adult; detection and
  reconstruction claims rest on separate results. Simulated IED-source recovery does not
  identify an epileptogenic zone or establish surgical benefit.

## 12. Head motion and OPM slippage (G4, bounded secondary extension; NEW) — `opmsquid.motion`, `scripts/g4_motion.py`
Static fit is studied first (sections 8-11). This extension adds a separate, time-varying analysis of
motion: a head-mounted array keeps its sensor-to-head geometry while its relationship to the residual
room field changes, and a cap can slip. This section bounds both mechanisms on the G3B arrays and
noise conventions (Neuromag at top contact; the refitted dense OPM array; intrinsic + brain noise;
10-nAm cortical-normal dipoles; area-weighted medians over cortical targets) for the adult and the
24- and 12-month templates (O'Reilly et al. 2021, built from the Neurodevelopmental MRI Database of
Richards et al. 2016). It does not establish motion robustness.
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
computed at ed852b2; medians over cortical targets of the adult and the 24- and 12-month templates;
dB of detectability, dense OPM or Neuromag combined, intrinsic + brain noise)
* A. Neuromag, head displaced in the fixed helmet, template of the reference position: 2 mm costs
  0.04-0.14 dB, 5 mm 0.28-0.56 dB, 10 mm down 1.64-1.72 dB (5-8 % of the cortex losing more than 3
  dB) and 10 mm sideways or forward/back 1.37-1.46 dB (templates only: at top contact the adult's
  head has no room); 5-deg rotations 0.22-0.84 dB, 10-deg rotations 0.91-2.31 dB (pitch and roll of
  the templates 1.90-2.31 dB, with 34-41 % of their cortex losing more than 3 dB). With the displaced
  geometry known (ideal movement compensation) the loss is at most 0.39 dB (10 mm away from the
  helmet top). 18 of the 81 displacements are infeasible at top contact (translations towards the
  helmet wall and some pitch and roll rotations; no upward displacement was tested). Dense OPM cap
  slipped, template of the reference geometry: 1 deg (median sensor shift 1.3-1.9 mm) costs
  0.01-0.06 dB, 3 deg (3.8-5.6 mm) 0.12-0.46 dB; with the slip known at most 0.07 dB. A 3-deg slip
  lifts 11-84 sensors (largest lift 7.5 mm, adult, about x; v2 up to 24 mm), most on the adult's less
  spherical head.
  Per mm of sensor-to-head displacement the uncompensated losses of the two systems are similar; the
  head-mounted array is unaffected by head displacement itself, which costs the fixed helmet up to
  1.7 dB at 10 mm and 2.3 dB at 10 deg without compensation.
* B. In-band rotation, artefact outside the noise model; thresholds in deg RMS per axis for a unit
  field (divide by the residual field in nT or nT/m), 1-dB loss on the median curve over 32 draws,
  ranges over the three anatomies (the 10th-90th percentiles of the per-draw thresholds, a draw
  that does not reach the level within the tested rotations counting as beyond them, lie 2-23 %
  below and 4-23 % above them): no correction 0.016-0.023 deg (uniform) and 0.14-0.21 deg (gradient); homogeneous
  projection: the uniform term is removed exactly with perfect calibration, 0.42-0.46 deg with 1-deg/
  1-% calibration errors and 0.14-0.15 deg with 3 deg/3 %, while the gradient term is not removed
  (0.14-0.21 deg, as without correction); 8-term projection: uniform 0.40-0.42 and 0.13-0.14 deg,
  gradient 3.2-3.9 and 1.1-1.3 deg (1 deg/1 % and 3 deg/3 %). The OPM falls to Neuromag's static
  detectability (D = 0) at 0.68-1.5 times these rotations, earliest for the adult, whose static D
  after the projections is smallest (+0.67 to +0.87 dB, against +1.25 to +1.84 dB for the
  templates). With the head origin instead of the neck as pivot the thresholds after the homogeneous
  projection change by less than 2 % (not at all with perfect calibration: the translation the neck
  pivot adds is uniform and removed),
  and after the 8-term projection they rise by 13-22 % (or beyond 5 deg): with calibration errors the
  projection leaks part of the translation term. With the artefact part of the noise
  model (oracle), the loss stays below 0.15 dB up to 5 deg in every case. Artefact per channel for 1
  deg RMS in the unit field: 15 pT without correction (uniform), 0.23-0.24 pT after the homogeneous
  projection with 1 deg/1 % errors, 23-28 fT for the gradient term after the 8-term projection. For
  example, with 1-deg/1-% calibration and the 8-term projection the dense array loses 1 dB at about
  0.4 deg RMS of in-band head rotation in a 1-nT residual field and 0.04 deg in 10 nT; without any
  correction at 0.02 deg in 1 nT.
* B'. Exact rigid motion (adult; 5-deg drift, 0.05-deg in-band jitter, 2 nT and 5 nT/m, 1 deg/1 %):
  per channel the exact in-band artefact is within 0.3 % of the linear prediction for 90 % of the
  channels (largest deviation 2.7 %, after the homogeneous projection; medians 1,395, 335 and 24 fT RMS after no, the homogeneous and
  the 8-term correction, against 89 fT of intrinsic OPM noise in the band); the drift moves each
  sensor's operating point by up to 278 pT (median 130 pT).
* Reading. A head-mounted array removes the geometry error that head motion causes in a fixed
  helmet, and converts motion into a field artefact whose size scales with the residual field and
  whose removal depends on calibration (or on modelling it from the data or from measured motion).
  These are bounds on two mechanisms in declared conditions, not an estimate of motion robustness
  in children.
