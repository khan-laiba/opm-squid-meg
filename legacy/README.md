# Legacy work (before the adult-to-pediatric study)

Scripts and outputs written on 2026-09-29, before `GOAL.md` defined the adult-to-pediatric
study. They are kept unchanged. The only change is that they were moved from the repository
root into this folder on 2026-09-30. Absolute paths stored inside their JSON outputs still
refer to the original root location.

Run them from this folder with the project environment, e.g.
`../.venv/bin/python replicate_figure3.py`. Outputs are written next to the scripts.

| Item | What it is | Status | Evidence |
|---|---|---|---|
| `replicate_figure3.py`, `figure3_output/` | Reproduction of Jas et al. (2026) Fig. 3: radial point magnetometers at 0 and 18 mm around a 95/80-mm sphere, 30 nAm tangential dipole, eta = 3 | **Reference reproduction, validated.** Reused by G1A. | `verify_figure3.py`: 7/7 checks pass against the published raster. MNE forward vs Biot-Savart within 6e-8 relative; interpolated peak vs Eq. 1 within 1e-7. Independent reviewer approved (2026-09-29). |
| `verify_figure3.py` | Pixel-level comparison of the replica with the figure embedded in the preprint PDF | Validation tool | Negative-control self-test passes. Needs the reference PDF in the repository root (not committed). |
| `simulate_extended_sources.py`, `extended_sources_output/` | Fig. 3 with spherical-cap patches (rho = 0-20 mm, 30 nAm or fixed density) instead of the dipole | **New study choice (not in the paper), validated.** Idealized sphere model. | Independent verifier reproduced every d_eq (27.67, 27.13, 25.55, 23.02, 19.73 mm) and peak to <= 6e-6. Robustness fixes approved (2026-09-30). |
| `simulate_realistic_head.py`, `realistic_head_output/` | fsaverage cortex, 1-layer BEM, conformal point-magnetometer arrays on the MRI scalp (OPM 0 mm, SQUID 18 mm), constant eta | **Idealized baseline, not a physical sensor system.** The 28,112 "sensors" per array are dense field-sampling points, not independent channels. SQUIDs are conformal at a constant 18 mm; the noise ratio eta is constant. | Independent verifier reproduced every number (peaks to 4.4e-6 with its own code); sampling and BEM-refinement convergence documented in the script docstring. |

## How these results may be used in the new study

- The Fig. 3 replica is the G1A reference reproduction; its recovered hidden parameters are
  documented in `replicate_figure3.py` and `docs/provenance_register.md`.
- The sphere-patch and conformal-array results are an **idealized baseline** only. They
  must not be presented as a comparison of physical OPM and Neuromag systems: real SQUID
  helmets are rigid (scalp gaps vary by region), real sensors have finite coils and
  intrinsic noise, planar gradiometers are not point magnetometers, and noise is not a
  constant ratio of two scalars.
- The conformal-array analysis showed that, at eta = 3, the equal-SNR depth statistic d50
  for fsaverage was within about 1 mm of the sphere model. That holds only at eta = 3 and
  only for idealized arrays; it is a legacy claim that the realistic-array study (G2)
  re-examines rather than assumes.
