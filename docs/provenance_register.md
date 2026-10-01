# Parameter, provenance and assumptions register

Every numerical choice used by the study, with its origin. Categories never mix:

- **J / HU / GO**: reported in Jas et al. 2026 / Hunold et al. 2016 / Goldenholz et al. 2009
  (page numbers refer to `docs/literature/*.md`).
- **HW**: hardware specification or vendor/site documentation.
- **R**: recovered by measurement from a published figure, not stated in the text.
- **A**: new study assumption (with rationale and the sensitivity analysis that covers it).
- **D**: study decision (design choice, not a physical parameter).
- **U**: unresolved ambiguity (kept open, not silently repaired).

Status: `verified` (checked against the source or reproduced numerically), `transcribed`
(copied from GOAL.md; primary source not available here), `assumed`, `open`.

## Jas et al. 2026 (analytical benchmark)

| ID | Parameter | Value | Source | Status | Used in |
|---|---|---|---|---|---|
| J-h | Adult head radius h | 95 mm | GOAL.md; paper Table 1 "Adult" | verified (legacy replica) | G1A |
| J-b | Adult brain radius b | 80 mm | GOAL.md; paper Table 1 | verified | G1A |
| J-Q | Tangential dipole moment | 30 nAm | paper | verified | G1A |
| J-xi | Scalp standoff OPM / SQUID | 0 / 18 mm | paper | verified | G1A |
| J-eta | Relative RMS noise ratio sigma_OPM / sigma_SQUID | 1-6 (3 in Fig. 3) | paper | verified | G1A |
| J-eq1 | Peak radial field of a tangential dipole (Eq. 1) | closed form, see `opmsquid.sphere.bmax_radial` | paper Eq. 1 | verified vs independent Sarvas maximum (< 1e-8) | G1A |
| J-deq | Equal-SNR depth at eta = 3 | 27.665 mm (paper: "approximately 28 mm") | Eq. 3, computed | verified | G1A |
| J-opm-standoff | OPM sensing-centre distance from the helmet inner surface | ~7 mm (2-mm helmet shell + half of a 10-mm cell) | paper p. 10 | verified (extraction) | G2 (A-OPM-STANDOFF) |
| J-opm-cell | OPM vapour-cell size | 10-mm cube, FieldLine Gen2, single axis, normal component | paper p. 10 | verified (extraction) | G2 (A-OPM-CELL) |
| J-table1 | Head models (h, b) | newborn (55, 48), 1 y (70, 62), 8 y (85, 73), adult (95, 80) mm | paper Table 1, p. 7 | verified (extraction) | G3 size benchmark |
| J-noise-lit | Literature intrinsic noise | SQUID 2-5 fT/sqrt(Hz), OPM 7-30 fT/sqrt(Hz) (Brookes et al. 2022) | paper p. 4 | verified (extraction) | A-OPM-NOISE sweep range |
| J-eta-meas | Measured noise ratio (not stated) | ~4.6 (Fig. 8B, sigma_OPM / sigma_SQUID at the closest runs) | derived from the figure | derived | context |
| J-toy | Toy experiment (Fig. 6) | target r_Q = 0.6 b; noise dipoles 0.4 b and 0.8 b; all 30 nAm +y; SNR = ratio of peak fields; xi 0-60 mm | pp. 9, 16-17 | reproduced (G1A) | G1A |

## Hunold et al. 2016 (G1B)

| ID | Parameter | Value as printed | Page | Status |
|---|---|---|---|---|
| HU-anat | Anatomy | 2 subjects, FreeSurfer white surface, ~300k nodes | 1148 | replaced by MNE sample subject (311,994 valid vertices; mean vertex area 0.648 mm^2 vs 0.654 implied) |
| HU-bem | Head model | 3 shells, 5120 triangles, 0.33 / 0.0042 / 0.33 S/m (brain/skull/scalp), Galerkin BEM | 1148, 1151 | conductivities used as printed; MNE linear collocation (U-HU-bem) |
| HU-coil | Coil integration | 4-point per sensor | 1151 | study coils 9014/9024 = MNE 'normal' rule (0.6-0.8 % vs accurate) |
| HU-depth | Depth | distance to the nearest scalp BEM node | 1148 | as printed (2,562-node BEM scalp). Mesh-sensitive: with the dense scalp instead, 41 % of the dipoles change 5-mm bin (mean +4 mm; review finding) |
| HU-orient | Orientation | angle to the normal of the nearest inner-skull BEM node, 0-90 deg | 1148 | as printed (folded with abs(cos)) |
| HU-bins | Bins | 20-60 mm x 0-90 deg, 5 mm x 10 deg | 1150 | as printed |
| HU-dipole | Focal strength | 600 nAm peak | 1151 | as printed (paper-specific, not physiological) |
| HU-patch | Patches | >20 mm^2 (first exceed), +/-10 deg orientation window, totals 612-678 nAm (median 622) | 1149-1151 | as printed; density = 622 nAm / median area (U-HU-density) |
| HU-bg | Background | random 10 % of nodes, EEG bands weighted 0.4-0.6, +/-10 nAm, 6 s at 1 kHz, one realization | 1148-1151 | as printed where stated (U-HU-bands); effective level unresolved (U-HU-bglevel) |
| HU-noise | Sensor noise | omitted | 1158 | as printed (reference); intrinsic noise only in the labelled extension |
| HU-snr | SNR | channel with max noise-free amplitude; 1 s pre-onset baseline; 2 mean abs(hilbert) | 1151 | as printed (U-HU-numerator) |
| HU-thr | Threshold | 2.5 (visual, one clinician) | 1151, 1159 | reported, not used as a detector threshold |
| HU-counts | Sources per bin | 3783 dipoles (Fig. 3c) | 1151-1152 | reproduced exactly by stratified sampling |

## Goldenholz et al. 2009 (G1C)

| ID | Parameter | Value as printed | Page | Status |
|---|---|---|---|---|
| GO-focal | Focal dipole | 10 nAm, every vertex, cortical normal | 1078 | as printed |
| GO-patch | Patches | geodesic radius 10 / 16 mm (3 / 8 cm^2), 50 pAm/mm^2 | 1078 | as printed; Dijkstra on the full mesh, centroids = oct-6 vertices (U-GO-patch) |
| GO-bem | Conductivities | 0.3 / 0.06 / 0.3 S/m (brain/skull/scalp) | 1078 | run as printed AND with 0.006 (probable intended MNE default; U-GO-skull) |
| GO-eq1 | SNR | 10 log10[(a^2/N) sum_k b_k^2/s_k^2] | 1079 | as printed; mag, grad and pooled reported separately (U-GO-pool) |
| GO-noise | Modelled noise | independent cortex-normal sources, ~7 mm grid, s_k^2 = s_s^2 (AA^T)_kk | 1079-1080 | as printed; 7-mm Poisson-disk grid, 1,736 usable sources (U-GO-grid) |
| GO-cal | Calibration | per type median(recorded/(AA^T)), channel-weighted mean; s_s = 1.6-1.9 nAm | 1080 | as printed (EEG term absent); s_s = 2.68 nAm on our grid, 1.76 nAm normalised to 4,000 sources |
| GO-rec | Recorded noise | 2 min spontaneous, 0.5-100 Hz, magnetometer SSP | 1078-1079 | ADAPT: pre-stimulus baselines (-200-0 ms) of all 320 events of the sample task recording (317 kept after edge exclusion), sample SSP; restricting to the auditory/visual events (IDs 1-4) changes s_s by 0.6 % (review) |

## Recovered from the published Fig. 3 raster (legacy replica)

| ID | Parameter | Value | Status |
|---|---|---|---|
| R-sigma | Absolute SQUID noise sigma_SQUID | B_max(0.8 b, h + 18 mm) = 0.35458 pT (SNR_SQUID = 1 at d = 31 mm); sigma_OPM = eta sigma_SQUID | recovered; the paper gives only eta. Absolute SNR values depend on it, ratios and d_eq do not. |
| R-grid | Curve sampling | r_Q = linspace(0, b, 101)[1:] (d = 15 ... 94.2 mm) | recovered |
| R-marker | Dotted d_eq line position | 27.530 mm (deepest point with SNR_OPM > SNR_SQUID on d = linspace(15, 95, 250)) vs exact root 27.665 mm | recovered; authors' exact rule not identifiable (U-J1) |

## Recovered from Hunold et al. 2016 Fig. 6 (G1B; `scripts/digitise_hunold_fig6.py`)

| ID | Parameter | Value | Status |
|---|---|---|---|
| R-HU-bars | Scale-bar brackets in the embedded 200-ppi raster (serif to serif) | EEG 24.99 px = 100 uV; MM 25.36 px = 5 pT; GM 25.33 px = "100 pT" (read as pT/m); time axis 130.7 px/s | recovered |
| R-HU-baselines | Baseline SD of the drawn line (per-column centroid, display time < 0.90 s), per channel | MM 0631 3.08 px (0.61 pT), 0711 3.15 px (0.62 pT), 0741 2.91 px (0.57 pT); GM 0413 2.64 px (10.4 pT/m), 0412 2.28 px (9.0 pT/m), 0423 2.66 px (10.5 pT/m) | recovered; compared with our traces drawn the same way (`hunold.rendered_centroid_sd`) |
| R-HU-spike | Spike-window extreme of the superficial tangential dipole | MM 0631 23.9 px = 4.7 pT (noisy trace) | recovered; context for U-HU-bglevel |

## Hardware

| ID | Parameter | Value | Source | Status | Used in |
|---|---|---|---|---|---|
| HW-mag-noise | Neuromag TRIUX magnetometer typical white noise | 3.5 fT/sqrt(Hz) = 3.5e-15 T/sqrt(Hz) | TRIUX specification image (user-supplied) | transcribed from GOAL.md; image not on disk (U-HW1) | G2 |
| HW-grad-noise | Neuromag TRIUX planar gradiometer typical white noise | 3.6 fT/(cm sqrt(Hz)) = 3.6e-13 T/(m sqrt(Hz)) | same | transcribed (U-HW1) | G2 |
| HW-18mm | Pickup coil to room-temperature (helmet) surface | 18 mm | same | transcribed. Not a scalp gap: measured scalp-to-coil distances for the sample subject are 25-41 mm (median 31 mm) | G2 |
| HW-T3 | MRN coil types | planar gradiometer 3014, magnetometer 3024 (T3) | mrn.org (read 2026-09-30) | verified | G2 |
| HW-MRN-floor | MRN shielded-room noise floor | "5-7 fT" (bandwidth unstated) | mrn.org | transcribed; not used as a spectral density | context |
| HW-geometry | Neuromag sensor geometry | MNE sample recording (MGH Vectorview, 306 channels; stored coils 3012 x 204, 3024 x 102) | MNE 1.13.2 sample dataset | verified (read) | G1B, G1C, G2 |
| HW-coildef | Coil integration | MNE 1.13.2 `coil_def.dat`, 'accurate' level (3014: 8 points, 26.39-mm coils, 16.80-mm baseline; 3024: 16 points, 21-mm coil) | MNE | verified | all |

## New study assumptions (A) and decisions (D)

| ID | Item | Value / rule | Rationale and coverage |
|---|---|---|---|
| D-T3 | Gradiometer coil variant | 3012 in the sample file replaced by 3014 | Same geometry and integration points in MNE 1.13.2 except the 1-point model (0.3 mm offset); matches MRN (HW-T3). |
| A-OPM-CELL | OPM sensing volume | 10-mm cube, custom coil 9901 (27-point Gauss at MNE 'accurate') | J-opm-cell; MNE forward = independent quadrature (< 1e-6). Effect at the field peak: -1.9 % (10 mm source distance), -0.4 % (15 mm), -0.04 % (30 mm). |
| A-OPM-STANDOFF | Sensing centre from helmet inner surface | 7 mm | J-opm-standoff |
| A-OPM-GAP | Helmet inner surface to MRI scalp | 0 mm baseline; swept 3 and 6 mm by moving the primary sites outward along their axes (v2; v1 rebuilt the arrays, 215/223/231 dense sites) | Hair and helmet fit are not modelled otherwise. |
| A-OPM-AXIS | OPM sensitive axis | normal of the smooth BEM head surface (5,120 triangles) averaged within 15 mm, single axis | The dense scalp mesh has non-manifold edges and some inward normals; its smoothed normals put 33 dense sites 30-132 deg off (G2 review). Now within 7.1 deg of the nearest head-surface vertex normal at every site (matched array <= 5.5 deg; one midline occipital dense site at 7.05 deg; 95 % of sites within 3.5 deg; `tests/test_g2.py`). Multi-axis OPMs are a secondary extension. |
| A-OPM-COVER | OPM coverage | scalp above the plane through LPA, RPA and a point 3 cm above the nasion (it slopes down towards the back, so the occipital scalp to ~3 cm below the fiducial plane is included); no site within 20 mm of the preauricular points. Dense arrays: scalp points > 2 mm inside or > 4 mm outside the smooth BEM head surface (ear canals and pinna folds; the pinna) and scalp points within 20 mm (MRI z) of the lowest point of the MRI head surface are not sites | Excludes face, ears and neck. The sample MRI is pitched ~36 deg relative to the head frame, so the flat cap where its head surface is cut at the edge of the field of view lies above the brow plane at the back; an earlier dense array had 8 sites on that cap (axes into the neck) and 2 on the pinna. Sample subject: 98 of 102 Neuromag sites pass the coverage rules and 97 remain after the clearance rules (A-OPM-CLEAR); the dense array's lowest scalp point is 38 mm below the matched array's (MRI z -105 vs -67 mm). |
| A-OPM-CLEAR | Physical and model clearance | sensing centre >= standoff - 1 mm (6 mm) from every MRI scalp point and >= 4 mm from the BEM head surface; every integration point of the 10-mm cell >= 1 mm outside the BEM head surface (v2; exact point-to-triangle distance with the inside test, the cell oriented as the forward model builds it in the head frame; the first v2 check built it in the MRI frame, with its in-plane axes rotated by up to ~100 deg); sites moved outward along their axis in 0.5-mm steps, at most 5 mm, else infeasible and dropped | Package collisions with the pinna, brow and occipital curvature, and no model boundary inside a cell (v1 had integration points up to 2.4 mm inside the head surface, where the 3-layer BEM field is unreliable; pre-freeze review F2). Sample: matched 97 of 102 sites (19 moved, 0.5-5 mm), dense 212 (50 moved); centres 5.3-11.5 mm (median 6.3 mm) above the BEM head surface; every integration point >= 1.000 mm outside it. |
| A-OPM-NOISE | OPM intrinsic noise | sweep 7-30 fT/sqrt(Hz) (white) | No single verified device specification (GOAL.md). |
| D-ANAT | Primary adult anatomy | MNE `sample` subject | Individual adult MRI with real helmet position; fsaverage secondary. |
| D-BADCH | Bad channels | channels marked bad in the sample recording (MEG 2443) are excluded wherever recorded or empty-room noise enters | Its pre-stimulus RMS is 23x the gradiometer median. |
| A-BEM-DIST | Usable cortical sources | white-surface vertices inside the inner skull and >= 4 mm from its 5,120-triangle mesh (91.3 % of the valid vertices, 7,661 of the oct-6 targets) | Refining the mesh to 20,480 triangles changes Neuromag gains by a median 14 % (90th pct 67 %, max 341 %) at 2-3 mm, 2.5 % (13 %, 142 %) at 3-4 mm, 0.8 % (3.6 %, 58 %) at 4-5 mm, 0.4 % (1.9 %) at 5-6 mm; OPM about half. Excluding the nearest-to-skull crowns makes the superficial OPM advantage conservative. |
| A-BEM-SKIN | 3-layer head surface | the 5,120-triangle BEM head surface subdivided once (20,480 triangles, same flat geometry) in every 3-layer model (v2) | On this head the field a few mm outside the coarse surface is not converged: refinement changed the field at OPM integration points by a median 0.7-0.8 % of the array field scale at 3-4 mm (95th pct 8-12 %) and the dense-array headline by -2.4 % on the v1 arrays, and by 0.66-0.74 % (95th pct 7.5-9.8 %) and -1.6 % on the v2 arrays (every integration point >= 1 mm outside; 1.147x coarse vs 1.129x refined); Neuromag gains 0.08 %. An exact 3-shell sphere test shows the BEM code is accurate near a regular surface (95th pct 0.45 % at 2 mm, 5,120 triangles). A second refinement is not feasible here (full BEM matrix); assuming the error falls as h^2 (mesh size), as on the sphere (an assumption, not verified on this head), the refined surface is within ~1 % for the headline. Scripts: study_opm_near_mesh.py, study_bem_sphere_accuracy.py, study_bem_skin_refinement.py. |
| A-G2-BAND | Common analysis band | 1-40 Hz, Butterworth order 4, zero phase, for the target, every noise term and the measured-noise calibration; ENBW of the composite response | Band sensitivity (`scripts/g2_band_sensitivity.py`): 1-10, 8-30, 30-80 Hz with the brain scale and room field recalibrated per band. |
| A-G2-BRAIN | Cortical background | independent, cortex-normal, area-scaled moment variance on a 7-mm Poisson-disk grid over usable vertices (1,755 sources); scale fitted once so the median good-gradiometer variance equals the measured (task baseline - empty room) level; correlated extension exp(-d/lambda), lambda 5 and 10 mm, re-fitted | Gradiometers are least affected by room and physiological artefacts; the magnetometer level is a validation, not a fit. |
| A-G2-ENV | Room field | 8-term external expansion (homogeneous + linear gradient) about (0, 0, 40) mm head, coefficient covariance from the empty-room recording; same room field for every array; head-position variants keep it in head coordinates | Residual higher-order and movement-related fields not modelled. |
| A-OPM-BW | OPM frequency response (band supplement) | first-order low-pass, 100-Hz corner, applied to the OPM target and brain/room terms (intrinsic noise at the output) | SERF OPMs of this class have ~100-150 Hz bandwidth; matters only in the 30-80 Hz band. |
| A-G2-OPMNOISE | Primary OPM noise | 15 fT/sqrt(Hz), the middle of the A-OPM-NOISE sweep, not a device value | The sweep is reported for every comparison. |
| A-G2-TARGET | Focal targets | 10-nAm dipoles, cortical normal, at the 7,661 usable oct-6 vertices (A-BEM-DIST) | Every SNR metric is linear in the moment. |
| A-G2-PATCH | Extended targets | geodesic radius 5/10/20 mm around each target, signed sum; fixed total 10 nAm (scalar moment) or fixed density 0.25 nAm/mm^2 | OPM/SQUID ratios do not depend on the convention; absolute values do. |
| A-G2-COVEST | Estimated covariance | Ledoit-Wolf estimate from 2 x 39 Hz x T independent samples, T = 10 and 60 s; plug-in matched filter evaluated under the true covariance | Practical counterpart of the oracle detectability. |
| A-G2-HEADPOS | Head-position variants | measured; +/-5 mm along each device axis; +/-5 deg pitch about the head origin; well fitted (up 6 mm, 20 mm from the nearest magnetometer); all at least 18 mm (Dewar spacing) from the coils | Source-blind (geometry only); OPM arrays are head-mounted and unchanged. |
| D-G2-COND | Noise conditions | intrinsic; intrinsic+brain; intrinsic+brain+env; projected (external 8-dim subspace removed by the same noise-weighted projection for every array; rank n - 8) | Headline conditions: intrinsic+brain and projected. |
| A-G4-TS | Time-domain noise (G4) | G2 model in the time domain, one realization shared by all arrays: independent 1/f background moments on the G2 grid (band variance = brain scale x area), room field synthesised with the empty-room cross-spectra (peaks near 12 and 30 Hz) and scaled to the fitted band covariance, white intrinsic noise; 1-40 Hz, decimated to 150 Hz | Spectral shape of brain noise (1/f) is a choice; the band covariance equals the G2 model. |
| A-G4-EVENT | IED-like events | Hunold spike-wave complex (digitised) through the analysis filter; stretch 0.75/1/1.5; focal dipoles and 10-mm patches (fixed total) at 10-320 nAm; 72 locations stratified by depth x orientation (never on the medial wall, 'unknown'; v2); one event every 2 s | Morphology and strength swept, not only visible spikes. |
| A-G4-DET | Detectors | whitener from 10 min of baseline null data; oracle (known topography, template and time; per-trial false-positive probability 0.001); practical scanner (template bank x 716-candidate cortical dictionary distinct from the targets; events = local maxima >= 0.25 s apart; thresholds for 1 and 0.2 false events per minute from 20 min of calibration null data, frozen, checked on 20 min held out); hit = statistic above threshold within +/-50 ms of the true peak | Event-level false alarms (per minute) are distinct from the oracle's per-trial false-positive probability. |
| A-G4-LOC | Localization (bounded) | 24 locations, focal and patch, 80 and 320 nAm; inverse with a 1-layer BEM and a 2-mm/2-deg coregistration error (8 random draws shared by all arrays; location i uses draw i mod 8, so every condition sees all 8 draws; v1 used event e mod 8, 2 draws per condition), noise covariance from 5 min of null data; MNE/dSPM (depth 0.8, SNR 3) on a 5-mm grid (3,821 usable vertices) excluding the true vertices; ECD with MNE fit_dipole (min_dist 5 mm); errors on the MRI through the analyst's (perturbed) transform; no GOF cut | Model and registration mismatch avoid an inverse crime; no SNR equalization across arrays. MNE's GOF is computed on whitened data, where noise adds about one unit per channel, so a fixed GOF threshold would favour arrays with fewer channels; GOF, chi2/dof and the 95 % confidence volume are descriptive. |
| D-G3-ANAT | Pediatric anatomy (G3B) | the 2-year infant template ANTS2-0Years3T (O'Reilly et al. 2021, NeuroImage 227:117682; built from the Neurodevelopmental MRI Database, Richards et al. 2016), `mne.datasets.fetch_infant_template('2yr')`, 391,756,700-byte archive from github.com/christian-oreilly/infant_template_paper (LGPL-2.1), native dimensions; the adult scaled by 85/95 and by the template/adult head circumference as size-only controls | Owner decision 2026-09-30 (download approved). No school-aged native anatomy was available (MNE packages infant templates up to 2 years; a school-aged dataset needs another download). The template is an average, not an individual; one pediatric anatomy, so no variability claim. Template files are not committed; figures derived from it cite both papers, and their redistribution (the Richards database has its own terms) is to be confirmed with the owner before any public release. |
| A-G3-SCALE | Size-only controls | every MRI-frame coordinate of the adult multiplied by the factor about the MRI origin (vertex-homologous; areas x factor^2); 85/95 = 0.895 (Jas Table 1 8-year/adult head radius) and the occipitofrontal-circumference ratio of the template to the adult | A scaled adult has adult proportions (tissue thickness, folding) at child size; comparing it with the template at the same head size separates size from anatomy. |
| A-G3-PLACE | Head placement in the fixed helmet | centred: the adult's measured device-to-head transform; top/back: moved along device +z / -y until the nearest magnetometer coil centre is 20 mm from the scalp; variants: +-5 mm along x, y, pitch +-10 deg, roll +-5 deg about the head origin, then top contact where there is room; feasible if every coil centre is >= 18 mm from the scalp. Primary: top | Source-blind (scalp and helmet only). 20 mm = HW-18mm + 2 mm, G2's 'well fitted' rule. Children are usually positioned against the top of the helmet; the centred and back placements and the variants bound the dependence. The adult's measured pose need not centre another head laterally: the template sits right of the helmet's midline (centred: left/right temporal gaps 48.9/33.3 mm); the 'x-centred' variant shifts each head along device x to equal left/right median gaps before the top contact, and 'top-18mm' is true contact at the Dewar spacing. |
| A-G3-COUNTERFACTUAL | Reduced-mismatch helmet | coil centres scaled by the child/adult head-circumference ratio about the head origin (centred placement); raised in 0.005 steps if a coil comes within 18 mm of the scalp | Mechanistic control holding coils and noise fixed; reproduces the adult's fit, scaled, around a scaled adult. Not a pediatric SQUID system. For the template the factor is limited by its lateral offset (right side), so the helmet is also scaled about the laterally centred head ('counterfactual_x-centred'). |
| A-G3-METRIC | Pediatric comparison | D = 20 log10(d_OPM / d_SQUID) [dB] per SQUID comparator; Delta = D_child - D_adult; area-weighted medians over homologous sources (vertex-wise for scaled controls; parcels and depth/orientation strata for the template); parcel bootstrap | GOAL G3: report D_child, D_adult and Delta together; a positive Delta is a relative gain under these matching assumptions, not by itself an OPM advantage in the child. |
| A-G3-USEFUL | Usefulness maps | known-topography detectability >= 5 at a reference moment (20, 50, 100, 200 nAm; maps at 100 nAm) | Operational choice, not a clinical standard; the moment needed for d = 5 is also reported by depth. |

## Unresolved ambiguities (U)

| ID | Item | Handling |
|---|---|---|
| U-J1 | Jas Fig. 3 d_eq marker at 27.53 mm vs Eq. 3 root 27.665 mm vs "28 mm" in the text | Both values reported; exact root used in analyses. The replica's 250-point grid rule is fitted to the raster (many grid sizes give 27.50-27.56 mm), not the authors' rule; recorded in `g1a_benchmark.json` (markers). |
| U-J2 | Fig. 4B d_eq printed 34 mm, drawn ~35.0 mm, exact 35.67 mm; Fig. 4C printed 19 mm, drawn ~19.8, exact 19.82 | Exact roots used; printed and drawn values listed in `results/g1a/g1a_benchmark.json`. Partly explained: the drawn markers are depths on the b/100 r_Q grid, and the printed values equal them truncated in SI floating point (0.034999... m -> 34, 0.019799... m -> 19); the same rule gives 26, not 28, for Fig. 3. |
| U-J3 | eta0 / eta1 printed 1.7 / 5.3; exact 1.6829 / 5.3086 (a deep crossing exists for 1.683 < eta < 1.7) | Exact values used. |
| U-J4 | Fig. 5B "normalized d_eq" = (d_eq - (h - b))/b (depth below the brain surface as a fraction of b), not d_eq/b; text speaks of brain *volume* | G3 reproduces the radial fraction as plotted and reports the true volume fraction separately (adult 40.4 %, newborn 87.2 % at eta = 3). |
| U-J5 | Fig. 5 uses s = h + 18 mm for every head (size-following shell) | Reproduced as the "size-following, constant-standoff benchmark"; the fixed-helmet experiment is separate (G3). |
| U-J6 | Fig. 6 depths: text 32/48/64 mm (these are r_Q = 0.4b/0.6b/0.8b) vs caption 63/47/31 mm; noise-dipole moments unstated (30 nAm reproduces the figure) | Caption values and 30 nAm used. |
| U-J7 | sigma_N20 defined as trial SD (p. 12) vs "standard error of mean" (Fig. 8) | Experimental context only (no recordings available). |
| U-J8 | Experimental OPM ~7 mm from the scalp (2-mm shell + half of a 10-mm FieldLine Gen2 cell) vs simulated OPM at xi = 0 | G1A keeps xi = 0 (paper model); G2 uses 7 mm (A-OPM-STANDOFF). |
| U-HW1 | TRIUX specification image not available on disk | Values transcribed in GOAL.md used and flagged; to verify when the image is supplied. |
| U-HU-trace | Hunold "dipole traces" (sulcal bottom to crown) cannot be replicated on new anatomy; algorithm partly unstated | Stratified sampling reproducing the paper's per-bin counts exactly (ADAPT). |
| U-HU-numerator | Spike amplitude in the SNR (noise-free vs noisy, peak vs peak-to-peak) unstated; the literal reading appears ~30-40 % below published values | Primary: noise-free peak-to-peak (Fig. 6 digitisation, one convention throughout: 0.87-1.10x the printed values; literal noisy peak 0.61-0.73x). Variants peak, noisy peak and noisy peak-to-peak are reported. Calibrated maps: noise-free p2p 0.78-0.91x (strong bins 0.87-0.92x, weak 0.73-0.90x); noisy p2p 0.97-1.04x (strong 0.91-0.97x, weak 1.01-1.16x); the paper lies between them, closer to the noisy p2p. The primary is not changed to fit the maps. |
| U-HU-bands | EEG band edges, weights, filter, normalisation unstated | 0.5-4 / 4-8 / 8-13 / 13-30 / 30-45 Hz, weights 0.6/0.55/0.5/0.45/0.4, Butterworth 4 zero-phase, unit-RMS per band, peak-normalised to 10 nAm per dipole. Noise generated 3 s longer on each side and cropped, so the 6 s are stationary as the paper states (without the padding, filter edge transients hold the maximum of ~90 % of the dipoles and the stationary part shrinks ~1.6x). |
| U-HU-bglevel | Effective background level at the sensors. With the text's +/-10 nAm per dipole (stationary), our expected rendered baselines at the Fig. 6 channels are 2.1x (MM) and 2.4x (GM) the paper's (R-HU-bars, R-HU-baselines; mean over 20 background realizations). The MM example spike agrees (6.5 pT) but our GM spike is 1.6-1.8x the paper's; normalisation, filtering, head position and anatomy cannot be told apart. | Both reported. `as_specified`: the text's level (bin means 0.37-0.43x the paper). `fig6_calibrated`: one scalar on all background moments, 0.47 (paper's MM baselines at 0631/0711/0741 over our expected baselines; one realization alone gives 0.44-0.52, 5-95 %, and an earlier single-realization run gave 0.43); other channel choices and statistics give 0.42-0.54, i.e. p2p bin means 0.68-1.01x. Not fitted to the maps. |
| U-HU-density | Patch density (implied ~30 nAm/mm^2; the discussion mentions 100 nA/mm^2) | Density set so that the median patch total is 622 nAm. |
| U-HU-onset | Spike onset within the 6-s epoch | 3.0 s; analytic signal of the whole trace, cropped to the 1-s baseline. |
| U-HU-bem | Galerkin vs MNE linear collocation BEM | Small differences expected (paper: MEG weakly sensitive to skull conductivity). |
| U-GO-skull | Printed skull 0.06 S/m (1:5) vs probable 0.006 (MNE default; manual v2.5 prints 0.06 beside "1/50") | Both run and labelled; MEG differences reported. |
| U-GO-pool | Whether Eq. 1 pooled all 306 MEG channels | Magnetometers (N=102), gradiometers (N=203) and pooled (N=305; MEG 2443 excluded, D-BADCH) reported separately. |
| U-GO-grid | Noise-source count and decimation unstated | 7-mm Poisson-disk grid on 3-D distance over usable vertices: 1,736 sources, sparser than the ~4,000 of a 7-mm MNE surface grid (A3 of the extraction); s_s is compared after scaling by sqrt(M / 4000). The review found the median pooled SNR within -0.19 to +0.05 dB for 15 of 17 alternative grids (M = 320-6,568); the other two contained a near-skull artefact source, now excluded by A-BEM-DIST. |
| U-GO-patch | Centroid set, surface and element weighting for patches unstated | oct-6 vertices as centroids; white surface; moment = density x vertex area; signed sum. |
| U-OPM-PACK | OPM package footprint (no verified device data) | Minimum sensing-centre spacing 17 mm (10-mm cell in a ~12-17 mm package). `opm_dense` is a greedy farthest-point array (212 sites on the sample head), not proven maximal (on an earlier candidate region the review found 2-4 % more sites with greedy insertion or random sequential adsorption). 306 single-axis channels appear infeasible (best spacing found 13.7-14.2 mm < 17 mm) but this is not proven. The channel-budget control `opm204` takes 204 dense sites by farthest-point sampling. The dense/matched contrast mixes channel count and coverage (the dense array reaches lower at the back than the Neuromag helmet; A-OPM-COVER). |
