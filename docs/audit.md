# G0 audit of existing work (2026-09-30)

Repository instructions and problem statement (recorded 2026-10-02): the repository holds no
instruction file of its own besides this study's specification and `PLAN.md` (no CLAUDE.md, AGENTS.md or
CONTRIBUTING file).

Scope: everything present before the adult-to-pediatric study, i.e. the scripts and outputs
now in `legacy/`, the anatomy and data on disk, and the sibling folder `/Volumes/T9/OPM_SQUID`
(a separate, parallel implementation, read but not modified).

## 1. Legacy scripts and outputs (`legacy/`)

| Component | Validated claims (evidence) | Assumptions and unverified items | Reuse in the study |
|---|---|---|---|
| `replicate_figure3.py` + `figure3_output/` | Reproduces Jas et al. Fig. 3 from MNE sphere forward simulations. MNE forward matches Biot-Savart to 5.9e-8 relative; the interpolated peak matches Eq. 1 to 9.9e-8. `verify_figure3.py` passes 7/7 checks against the published raster (curve median error 0.08-0.12 px, identical d_eq line, 156/156 text and arrow elements matched). Independent code review approved. | Parameters **recovered from the raster, not stated in the paper**: sigma_SQUID = B_max(0.8 b, h + 18 mm) = 0.35458 pT; curves sampled at r_Q = linspace(0, b, 101)[1:]; the dotted d_eq marker at 27.530 mm (the deepest point of a 250-point grid with OPM > SQUID), against the exact root of Eq. 3 at 27.665 mm. The authors' exact marker rule cannot be identified uniquely. Fonts and artwork offsets are cosmetic. | G1A reference reproduction of Fig. 3; `opmsquid.sphere` re-implements the physics independently and is tested against it. |
| `verify_figure3.py` | Negative-control self-test passes. | Needs the reference PDF, which is not committed. | G1A validation tool. |
| `simulate_extended_sources.py` + `extended_sources_output/` | Spherical-cap patches (rho 0-20 mm): an independent verifier reproduced d_eq = 27.67/27.13/25.55/23.02/19.73 mm and every peak to <= 6e-6. Robustness fixes re-verified. | **New study choice, not in the paper.** Caps carry a uniform +y current (tangential at the centre) with the whole cap at one depth; this is an idealization of cortex. | Sphere-model reference curves for patch-size comparisons (idealized). |
| `simulate_realistic_head.py` + `realistic_head_output/` | fsaverage (ico-5, 20,484 sources, cortical-normal), 1-layer BEM; the verifier reproduced every number. Sensor sampling converged (<= 0.3 % peak change); BEM refinement changes d50 by <= 0.03 mm. | **Idealized baseline, not a physical system.** Its 28,112 "sensors" per array are conformal field-sampling points; the SQUID array follows the scalp at a constant 18 mm; eta is constant; noise is not modelled explicitly; the helmet coverage plane is an assumption. Claims that hold only under these conditions: d50 ~ sphere d_eq within ~1 mm at eta = 3 (not at other eta), 22.5 % of cortex OPM-favoured at eta = 3. A claim that folded patches are "more OPM-favoured" was retracted after area-matched comparison. | Idealized baseline for G2/G3 comparisons; its cortical rendering, geodesic patch operator and depth definitions are reusable ideas (re-implemented with tests). |

No legacy configuration files or unit tests existed (CLI options only). No forward-model cache
existed; `realistic_head_results.npz` stores computed peaks, not forward solutions.

## 2. Anatomy and data on disk

| Item | Location | Notes |
|---|---|---|
| fsaverage (surfaces, 1- and 3-layer BEM surfaces, head, fiducials, ico-5 source space) | `/Volumes/T9/OPM_SQUID/data/subjects/fsaverage` (sibling folder; MNE config `SUBJECTS_DIR`) | Read-only use. Template anatomy: averaged, smoothed; not an individual. |
| MNE sample data (MNE 1.13.2 dataset) | `data/external/MNE-sample-data` (not committed) | Individual adult MRI subject `sample` (FreeSurfer surfaces, 3-layer BEM surfaces, oct-6 source space, head surface); MGH Neuromag Vectorview recording with real `dev_head_t`; head-MRI transform; empty-room recording and noise covariances. Coil types in the file: 204 x 3012 (planar gradiometer T1), 102 x 3024 (magnetometer T3). |
| Pediatric anatomy | none | Needed for G3 (see PLAN.md, Inputs). |
| Hunold/Goldenholz participant data, Jas SEF recordings | none | Reproductions are therefore adaptations. |

## 3. Sibling project `/Volumes/T9/OPM_SQUID` (context only)

A separate implementation of the same goal file, created independently. It contains a Fig. 3
reproduction, a spherical extension, and an anatomical extension (fsaverage, 3-layer BEM with
skull 0.006 S/m, 512 matched sensor sites, fixed-total and fixed-density normalizations) with
unit tests. Its results have **not** been verified here; they are not used as evidence in this
study. Its README states the same 27.665 mm equal-SNR depth.

## 4. Consequences for the plan

1. G1A can build on the validated Fig. 3 replica; Fig. 4, the absolute-noise documentation,
   the toy target/background experiment and the text/caption inconsistencies remain to do.
2. Everything physical in G2 (finite coils, rigid helmet, intrinsic and environmental noise,
   gradiometers) is new: no legacy result may be cited as a physical OPM-vs-Neuromag finding.
3. The legacy recovered hidden parameters stay labelled as recovered assumptions in the
   provenance register.
