# Scripts

The analysis, figure and fact scripts of the paper. Run them from the repository root with the project's
environment, for example `.venv/bin/python scripts/g2_adult_comparison.py`; the scripts that use the package add
`src/` to the import path themselves. `scripts/run_all.sh` runs them in the order the analyses depend on one another
(see the main `README.md`, "Reproducing the paper"). The paper itself is built by the scripts in `paper/`
(`export_figures.py`, `build_paper.py`, `build_html.py`; main `README.md`, "Reproducing the paper").

External data are read from `data/external/` and computed caches are written to `cache/` (neither is committed; the
locations can be changed with the environment variables `OPMSQUID_DATA` and `OPMSQUID_CACHE`). "MNE sample data"
below means `data/external/MNE-sample-data`, the anatomies are those of the MNE sample subject, the infant templates
(`data/external/infant_subjects/`) and the school-aged children (`data/external/school_subjects/`). The analyses
write their outputs to `results/`, and every JSON or CSV file there records the commit of the code that produced it.

`<head>` stands for one of the nine head models: `adult` (MNE sample subject), `school` and `size2yr` (the adult
scaled to school-age size and to the head circumference of the 24-month template), `infant2yr`, `infant18mo` and
`infant12mo` (the 24-, 18- and 12-month infant templates) and `childA`, `childB` and `childC` (the school-aged
children of OpenNeuro ds005234). Paper references: section, figure and table numbers of the main text; S1–S10 are the
sections, Figs. S1–S18 and Tables S1–S19 the figures and tables of the supplementary material.

## Figures and tables of the paper

Figures 8–10 show the adult and the four principal smaller heads; `paper/export_figures.py` draws them with the same
code restricted to these heads, and the versions with all nine heads are Figs. S16–S18. Tables are typeset from the
facts of `scripts/report_facts_*.py`; `paper/build/numbers_used.tsv` gives the source of every value. `<head>` is one of
the nine head models (`adult`, `school`, `size2yr`, `infant2yr`, `infant18mo`, `infant12mo`, `childA`, `childB`,
`childC`).

| Item | Content | Computed by | Result files | Drawn by (image in `results/report/`) |
|---|---|---|---|---|
| Fig. 1 | Sensor arrays on the adult head | `g2_adult_comparison.py`, exported by `export_g2_arrays.py` | `results/g2/g2_arrays.json` | `report_figures_clean.py` (`Figure_R15_arrays.png`) |
| Fig. 2 | The fixed and the fitted helmet | `g3b_pediatric_helmet.py`, `study_g3b_constant_gap.py`, exported by `export_g3b_geometry.py` | `results/g3b/g3b_geometry_sections.json` | `report_figures_clean.py` (`Figure_R12_geometry.png`) |
| Fig. 3 | Determinants of the adult comparison | `g2_adult_comparison.py` | `results/g2/g2_summary.json`, `results/g2/g2_targets.csv` | `report_figures_adult.py` (`Figure_R3_conditions.png`) |
| Fig. 4 | The noise model against measured Neuromag noise | `study_covariance_validation.py`, `study_noise_sensitivity.py` | `results/g2_covariance_validation/covariance_validation.json`, `results/g2_noise_sensitivity/noise_sensitivity_summary.json` | `report_figures_noise.py` (`Figure_R17_noise_checks.png`) |
| Fig. 5 | The adult comparison by depth | `g2_adult_comparison.py`, `study_g2_depth_bins.py` | `results/g2/g2_summary.json`, `results/g2/g2_depth_bins.json` | `report_figures_adult.py` (`Figure_R1_adult_depth.png`) |
| Fig. 6 | The adult comparison on the cortex | `g2_adult_comparison.py` | `results/g2/g2_targets.csv` | `report_figures_clean.py` (`Figure_R11_maps_adult.png`) |
| Fig. 7 | The adult comparison by brain region | `g2_adult_comparison.py` | `results/g2/g2_summary.json`, `results/g2/g2_targets.csv` | `report_figures_adult.py` (`Figure_R4_regions_adult.png`) |
| Fig. 8 | Smaller heads in the fixed adult helmet | `g3b_pediatric_helmet.py` | `results/g3b/g3b_summary.json`, `results/g3b/g3b_targets_<head>.csv` | `report_figures_pediatric.py` (A, B), `report_figures_clean.py` (C) |
| Fig. 9 | Helmet fit versus head size | `study_g3b_constant_gap.py`, `g3b_pediatric_helmet.py` | `results/g3b_constant_gap/g3b_constant_gap_summary.json`, `results/g3b/g3b_summary.json` | `paper/export_figures.py` |
| Fig. 10 | Simulated interictal spikes, pre-specified run | `g4_confirmatory.py`; exploratory ratios from `g4_epilepsy_pediatric.py --compare` | `results/g4_confirm/g4_confirm_summary.json`, `results/g4_confirm/g4c_<head>_summary.json`, `results/g4/g4_pediatric_comparison.json` | `report_figures_confirm.py` (`Figure_R16_confirm.png`) |
| Table 1 | Head models | `g3b_pediatric_helmet.py`, `study_children_qc.py` | `results/g3b/g3b_summary.json`, `results/g3b_templates_qc/children_qc.json` | facts |
| Table 2 | Main estimands | definitions (Section 2) | `results/g2/g2_summary.json`, `results/g3b/g3b_summary.json` | facts |
| Table 3 | The adult comparison by condition, comparator, measure, model variant and scenario | `g2_adult_comparison.py`, `g2_band_sensitivity.py`, `study_covariance_validation.py`, `study_noise_sensitivity.py` | `results/g2/`, `results/g2_covariance_validation/`, `results/g2_noise_sensitivity/` | facts |
| Table 4 | The principal smaller heads against the adult | `g3b_pediatric_helmet.py`, `study_g3b_constant_gap.py` | `results/g3b/g3b_summary.json`, `results/g3b_constant_gap/g3b_constant_gap_summary.json` | facts |
| Table 5 | Detection of simulated superficial spikes | `g4_confirmatory.py`, `g4_epilepsy_adult.py`, `g4_epilepsy_pediatric.py` | `results/g4_confirm/`, `results/g4/` | facts |

| Supplement | Content | Main scripts | Main result files |
|---|---|---|---|
| S1 | Estimands and conventions | as for the main text | `results/g2/`, `results/g3b/` |
| S2 | The spherical benchmark (Fig. S1) | `g1a_jas_benchmark.py` with `legacy/replicate_figure3.py`, `g3a_jas_size_benchmark.py` | `results/g1a/`, `results/g3a/` |
| S3 | The noise model, its validation and sensitivity analyses (Figs. S2–S5) | `g2_adult_comparison.py`, `g1b_hunold.py`, `g1c_goldenholz.py`, `study_covariance_validation.py`, `study_noise_sensitivity.py` | `results/g2/`, `results/g1b/`, `results/g1c/`, `results/g2_covariance_validation/`, `results/g2_noise_sensitivity/` |
| S4 | Head models and the MRI check (Fig. S6) | `prepare_school_subjects.py`, `study_school_anatomy.py`, `study_child_bem.py`, `study_children_qc.py`, `study_cortex_exclusion.py` | `results/g3b/`, `results/g3b_children_qc/`, `results/g3b_templates_qc/` |
| S5 | Helmet constructions (Figs. S7, S8) | `g3b_pediatric_helmet.py`, `study_g3b_constant_gap.py` | `results/g3b/`, `results/g3b_constant_gap/` |
| S6 | Regions and the noise floor across heads (Figs. S9–S11) | `g3b_pediatric_helmet.py`, `study_g3b_constant_gap.py` | `results/g3b/`, `results/g3b_constant_gap/` |
| S7 | The spike study (Figs. S12–S15) | `g4_epilepsy_adult.py`, `g4_localization.py`, `g4_epilepsy_pediatric.py`, `study_g4_matched_rate.py`, `g4_motion.py`, `g4_confirmatory.py`, `g4_confirm_exact_p.py` | `results/g4/`, `results/g4_confirm/` |
| S8 | Parameters and analysis choices | none beyond those above | `configs/`, `docs/methods.md`, `docs/provenance_register.md` |
| S9 | Clinical OPM studies of epilepsy and model analogs (Table S19) | none beyond those above | `docs/literature/epilepsy_opm_studies.json` |
| S10 | The 18-month template and the school-aged children (Figs. S16–S18) | as for Figs. 8–10 | `results/g3b/`, `results/g3b_constant_gap/`, `results/g4_confirm/` |

## Pipeline

| Script | What it does | Main inputs | Outputs | Paper |
|---|---|---|---|---|
| `run_all.sh` | Runs the unit tests, then the analysis, export and figure scripts below in a fixed order; `RESUME=1 scripts/run_all.sh` skips the steps already completed at the current commit | external data in `data/external/` | `results/` | the whole paper |

## Data preparation

| Script | What it does | Main inputs | Outputs | Paper |
|---|---|---|---|---|
| `fetch_school_subjects.py` | Downloads the files of the three school-aged children of OpenNeuro ds005234 (version 2.2.0), each from its exact S3 object version, and keeps a file only if its size and SHA-256 match the manifest; `--no-qc` leaves out the MRI-check files, `--check` checks the local files and downloads nothing | `configs/school_subjects_manifest.json`, `configs/school_subjects_qc_manifest.json` | `data/external/school_subjects/` | Section 2.2; S4, S10 |
| `prepare_school_subjects.py` | Prepares each child's anatomy: a three-layer BEM with a modeled skull, the dense MRI scalp, an oct-6 source space and fiducials transferred from the adult | `data/external/school_subjects/`, `configs/g3b_pediatric.toml` | `data/external/school_subjects/<subject>/bem/`, `results/g3b/school_subjects_preparation.json` | Section 2.2; S4.1 |
| `compute_fullres_forwards.py` | Lead fields of the adult at every usable white-surface vertex for seven array and BEM combinations; resumable, each matrix stored with a fingerprint of its array and BEM | MNE sample data | `cache/fullres/<job>.npy` | the Hunold and Goldenholz adaptations, the adult comparison, the spike study and the near-skull sensitivity analysis |

## Benchmarks

| Script | What it does | Main inputs | Outputs | Paper |
|---|---|---|---|---|
| `g1a_jas_benchmark.py` | The spherical-head analysis of Jas et al. (2026): Fig. 3 (regenerated by `legacy/replicate_figure3.py`), Fig. 4 and the Fig. 6 depth experiment, with independent numerical checks of the peak field | `legacy/replicate_figure3.py` | `results/g1a/` (`g1a_benchmark.json`, `g1a_curves.csv`, `fig3/`) | S2; Fig. S1 |
| `g3a_jas_size_benchmark.py` | The head-size benchmark of Jas et al. (2026) (their Table 1 and Fig. 5) in four concentric-sphere heads, and the same heads inside one fixed adult shell | none | `results/g3a/` (`g3a_size_benchmark.json`, `g3a_deq.csv`) | S2 |
| `g1b_hunold.py` | Adaptation of the depth–orientation spike simulations of Hunold et al. (2016), MEG part, to the adult, with a site-matched OPM array | `configs/hunold_reference.toml`, MNE sample data, `cache/fullres/` | `results/g1b/` (`g1b_summary.json`, `g1b_sources.csv`) | S3.2 |
| `digitise_hunold_fig6.py` | Digitizes Fig. 6 of Hunold et al. (2016) for the absolute calibration of the background; the constants in `src/opmsquid/hunold.py` were copied from its output | the published PDF (not redistributed) | printed values only | S3.2 |
| `g1c_goldenholz.py` | Adaptation of the cortical SNR maps of Goldenholz et al. (2009), MEG part, with an OPM extension | `configs/goldenholz_reference.toml`, MNE sample data, `cache/fullres/` | `results/g1c/` (`g1c_summary.json`, `g1c_oct6_values.csv`) | S3.2 |

## Adult comparison

| Script | What it does | Main inputs | Outputs | Paper |
|---|---|---|---|---|
| `g2_adult_comparison.py` | The adult comparison: Neuromag (306 channels) against the site-matched and the dense OPM array; known-topography detectability with sensor, brain and room noise, peak-channel and mean-power SNR, plug-in covariances, sensitivity and convergence checks | `configs/g2_adult.toml`, MNE sample data, `cache/fullres/` | `results/g2/g2_summary.json`, `g2_targets.csv`, `g2_patch_targets.csv` | Sections 3.1, 3.3; Figs. 3, 5–7; Tables 2, 3; S1, S3 |
| `study_bem_sphere_accuracy.py` | Exact test of the three-layer BEM field close to the outer surface, in concentric spheres meshed like the adult's BEM | none | `results/g2/bem_sphere_check.json` | Section 2.2; S8 |
| `study_opm_near_mesh.py` | Convergence of the BEM field at the integration points of the OPM cells (head surface of 5,120 against 20,480 triangles) | MNE sample data | `results/g2/near_mesh_check.json` | S8 |
| `study_bem_skin_refinement.py` | The adult ratios with the head surface refined and with a single-compartment BEM | MNE sample data | `results/g2/bem_skin_refinement.json` | Section 2.2; S8 |
| `study_head_surface_effect.py` | The adult ratios with the BEM head surface on the MRI scalp against the stored outer skin, separating the conductor, the OPM sites and the sensitive axes | MNE sample data | `results/g2/head_surface_effect.json` | Section 2.3; S3, S8 |
| `g2_report.py` | Text summary of the adult comparison per sensor configuration, from the stored outputs | `results/g2/g2_summary.json`, `results/g2/g2_band_sensitivity.json` | `results/g2/G2_report.md` | not used by the paper build |
| `study_g2_depth_bins.py` | 95 % parcel-bootstrap interval of the adult ratio in every 5-mm depth bin, from the stored per-target values | `results/g2/g2_targets.csv`, `g2_targets_depth.csv`, `g2_summary.json` | `results/g2/g2_depth_bins.json` | Section 3.3; Fig. 5 |
| `export_g2_arrays.py` | The adult's sensor arrays and head surface, rebuilt by the comparison's own code and checked against its summary, for drawing | MNE sample data, `results/g2/g2_summary.json`, `results/g3b/g3b_geometry_sections.json` | `results/g2/g2_arrays.json` | Fig. 1 |
| `export_g2_target_metrics.py` | Full-precision detectability of every adult target, with explicit flags where an OPM array is ahead | MNE sample data, `results/g2/g2_targets.csv`, `g2_summary.json` | `results/g2/g2_targets_metrics.csv` | Section 3.1 (target counts) |

## Noise-model checks and sensitivity analyses

| Script | What it does | Main inputs | Outputs | Paper |
|---|---|---|---|---|
| `g2_band_sensitivity.py` | The noise model rebuilt in other frequency bands, with the OPM sensor response, and a simulation check of the rule that converts spectra to variances | `configs/g2_adult.toml`, MNE sample data | `results/g2/g2_band_sensitivity.json` | Table 3; S3; Fig. S2 |
| `study_covariance_validation.py` | The noise model against the measured noise covariance of the Neuromag sample recording: channel variances, spatial correlation, eigenstructure, sub-bands, Neuromag's detectability with the measured covariance (finite-sample bias removed with surrogates) and scenarios A–E for the OPM | `configs/g2_adult.toml`, MNE sample data | `results/g2_covariance_validation/` (`covariance_validation.json`, `covariance_validation_targets.csv`) | Sections 2.4, 3.2; Fig. 4; Table 3; S3.4, S3.5; Figs. S3, S4 |
| `study_noise_sensitivity.py` | What the noise model leaves out: near-skull cortex, colored OPM noise, cardiac and ocular sources; the OPM white-noise sweep with the break-even levels; joint runs | `configs/g2_noise_sensitivity.toml`, `configs/g2_adult.toml`, MNE sample data, `cache/fullres/` | `results/g2_noise_sensitivity/noise_sensitivity_summary.json` | Sections 2.4, 3.2; Fig. 4; Table 3; S3.6; Fig. S5 |

## Smaller heads and helmets

| Script | What it does | Main inputs | Outputs | Paper |
|---|---|---|---|---|
| `study_school_anatomy.py` | Checks behind the children's anatomy: depth of the inner skull below the scalp, the fiducial transfer tested on the templates, the talairach transforms, the thinnest layers | anatomies in `data/external/`, `configs/school_subjects_manifest.json` | `results/g3b/school_anatomy_checks.json` | Section 2.2; S4, S10 |
| `study_child_bem.py` | The children's modeled skull applied to the adult, whose segmented skull is known, and its effect on the adult ratios | MNE sample data | `results/g3b/child_bem_validation.json` | S4 |
| `g3b_pediatric_helmet.py` | The adult and the eight smaller heads in the fixed adult helmet at source-blind placements (top contact primary) and in a helmet scaled with the head, with the OPM arrays refitted to each head; D and Δ with parcel-bootstrap intervals | `configs/g3b_pediatric.toml`, `configs/g2_adult.toml`, anatomies in `data/external/` | `results/g3b/g3b_summary.json`, `g3b_targets_<head>.csv`, `G3B_report.md` | Sections 2.2, 2.3, 3.4; Fig. 8; Tables 1, 4; S5, S6, S10; Figs. S7–S11, S16 |
| `study_g3b_constant_gap.py` | Each head in a helmet fitted at the adult's gap (and at the adult's top-contact gap): Δ, within-head contrasts, the helmet-fit interaction, placement bands and regions | `configs/g3b_pediatric.toml`, `configs/g2_adult.toml`, `results/g3b/g3b_summary.json` | `results/g3b_constant_gap/` (`g3b_constant_gap_summary.json`, `targets_<head>.csv`) | Sections 2.3, 3.5; Figs. 2, 9; Table 4; S5, S6, S10; Fig. S17 |
| `study_children_qc.py` | MRI quality check of the children's surfaces (with `--anatomies`, of other heads, such as the 18- and 12-month templates) against their T1-weighted images | T1 volumes and head masks in `data/external/`, `configs/school_subjects_qc_manifest.json` | `results/g3b_children_qc/`, `results/g3b_templates_qc/` (`children_qc.json`, figures) | Section 2.2; Table 1; S4.3, S10 |
| `study_cortex_exclusion.py` | The share of cortex that the 4-mm source rule leaves out in each head, overall and near the scalp | anatomies in `data/external/` | `results/g3b/g3b_cortex_exclusion.json` | Section 2.2; S4, S10 |
| `export_g3b_geometry.py` | The helmet geometry for drawing (scalp and white-surface sections, coil centers of the fixed, scaled and fitted helmets, dense OPM sites), exported from the stored state of the smaller-head run and checked against both summaries | `cache/g3b/state.pkl`, `results/g3b/g3b_summary.json`, `results/g3b_constant_gap/g3b_constant_gap_summary.json`, anatomies in `data/external/` | `results/g3b/g3b_geometry_sections.json` | Fig. 2 |

## Spikes

| Script | What it does | Main inputs | Outputs | Paper |
|---|---|---|---|---|
| `g4_epilepsy_adult.py` | The exploratory spike study in the adult: time-domain simulation with the adult noise model, oracle and practical scanning detectors, detection curves and S50 by depth band | `configs/g4_epilepsy.toml`, `configs/g2_adult.toml`, MNE sample data, `cache/fullres/` | `results/g4/g4_adult_summary.json`, `g4_adult_events.csv` | Section 2.5; Table 5; S7.1, S7.2; Fig. S13 |
| `g4_localization.py` | Bounded localization of the simulated spikes: dipole fits and dSPM through a single-compartment BEM with a 2-mm, 2° coregistration error | `configs/g4_epilepsy.toml`, MNE sample data | `results/g4/g4_localization_summary.json`, `g4_localization_events.csv` | Sections 2.5, 3.6; S7.3; Figs. S14, S15 |
| `study_g4_vs_g2.py` | Spike detection set against the detectability at the same locations | `results/g2/g2_targets.csv`, `results/g4/g4_adult_summary.json` | `results/g4/g4_g2_consistency.json` | S7.2 |
| `g4_epilepsy_pediatric.py` | The exploratory detection and localization studies on a smaller head at top contact (`<head> --detection --localization`), and the comparison with the adult (`--compare`) | `configs/g4_epilepsy.toml`, anatomies in `data/external/` | `results/g4/g4_<head>_summary.json`, `g4_<head>_events.csv`, `g4_localization_<head>_*`, `g4_pediatric_comparison.json`, `G4_pediatric_report.md` | Fig. 10; Table 5; S7, S10; Fig. S12 |
| `study_g4_matched_rate.py` | The exploratory detection re-evaluated with every threshold matched to 1 false event per minute on the held-out null data | `cache/g4/<head>_state.pkl` | `results/g4/g4_matched_rate.json`, `G4_matched_rate_report.md` | S7.2 |
| `study_g4_fit_failures.py` | Localization failures under declared criteria (gross errors, unconstrained dipole fits), counted from the per-event tables | `results/g4/g4_localization*_events.csv` | `results/g4/g4_fit_failures.json`, `G4_fit_failures_report.md` | S7.3 |
| `g4_motion.py` | Head motion and OPM slippage: sustained displacements in the fixed helmet, and in-band motion of the OPM array through a residual room field | `configs/g4_motion.toml`, `configs/g3b_pediatric.toml`, anatomies in `data/external/` | `results/g4/g4_motion_summary.json`, `G4_motion_report.md` | S7.4 |
| `g4_confirmatory.py` | The pre-specified spike run: new locations, seeds and null data, five noise replicates and a declared detector-mismatch variant; sign-flip tests Holm-adjusted over the anatomies (`<head> --no-combine` per head, then `--check-endpoint-code --combine-only`) | `configs/g4_confirmatory.toml`, `configs/g4_epilepsy.toml`, anatomies in `data/external/` | `results/g4_confirm/g4c_<head>_summary.json`, `g4c_<head>_events.csv`, `g4_confirm_summary.json`, `g4c_endpoint_code_check.json` | Sections 2.5, 3.6; Fig. 10; Table 5; S7.5, S10; Fig. S18 |
| `g4_confirm_exact_p.py` | Exact sign-flip p values of every comparison of the pre-specified run, with their Holm adjustment, from the stored per-location differences | `results/g4_confirm/` | `results/g4_confirm/g4_confirm_exact_p.json` | Section 2.5; S7.5, S10 |

## Figures and facts for the paper

The figure scripts draw from stored results only. `paper/export_figures.py` calls them again to export the print
figures, and every number of the manuscript is a fact of `report_facts.py` (main `README.md`, "Reproducing the paper").

| Script | What it does | Main inputs | Outputs | Paper |
|---|---|---|---|---|
| `report_style.py` | Shared look of the figures: open fonts, one color and marker per array and per head (Okabe–Ito palette) | none | none (imported) | all figures |
| `report_figures_adult.py` | Figures of the adult comparison | `results/g2/g2_summary.json`, `g2_band_sensitivity.json`, `g2_targets.csv`, `g2_depth_bins.json` | `results/report/Figure_R1_adult_depth.png`, `Figure_R2_noise_model.png`, `Figure_R3_conditions.png`, `Figure_R4_regions_adult.png`, `figures_adult.json` | Figs. 3, 5, 7, S2 |
| `report_figures_clean.py` | The spherical benchmark, the cortical maps, the helmet geometry and the adult's arrays (`--figures` draws a subset) | `results/g1a/`, `results/g2/` (`g2_summary.json`, `g2_targets.csv`, `g2_arrays.json`), `results/g3b/` (`g3b_targets_<head>.csv`, `g3b_geometry_sections.json`), `results/g3b_constant_gap/g3b_constant_gap_summary.json`; for the maps, inflated surfaces from `data/external/` | `results/report/Figure_R0_sphere.png`, `Figure_R11_maps_adult.png`, `Figure_R12_geometry.png`, `Figure_R13_maps_heads.png`, `Figure_R14_maps_scaled.png`, `Figure_R15_arrays.png`, `figures_clean.json` | Figs. 1, 2, 6, 8C, S1, S10, S16 |
| `report_figures_noise.py` | The checks of the noise model | `results/g2_covariance_validation/covariance_validation.json`, `results/g2_noise_sensitivity/noise_sensitivity_summary.json`, `results/g2/g2_summary.json`, `docs/methods.md` | `results/report/Figure_R17_noise_checks.png`, `figures_noise.json` | Fig. 4 |
| `report_figures_pediatric.py` | The smaller heads and the exploratory spike run | `results/g3b/g3b_summary.json`, `g3b_targets_<head>.csv`, `results/g4/g4_<head>_summary.json` | `results/report/Figure_R5_regions_heads.png` to `Figure_R10_spikes.png`, `figures_pediatric.json` | Figs. 8A–B, S7–S9, S11, S12, S16 |
| `report_figures_qc.py` | The figure of the MRI quality check, recomputed with the check's own functions and compared with its stored records | `results/g3b_children_qc/children_qc.json`, `configs/g3b_pediatric.toml`, T1 volumes in `data/external/` | `results/report/Figure_S_children_qc.png`, `figures_qc.json` | Fig. S6 |
| `report_figures_supplement.py` | Localization effect sizes, joint detection and localization, and the adult's detection curves | `results/g4/g4_localization*_summary.json`, `results/g4/g4_localization*_events.csv`, `results/g4/g4_adult_summary.json`, `results/g2/g2_summary.json`, `configs/g4_epilepsy.toml` | `results/report/Figure_S_localization_effects.png`, `Figure_S_joint_detection_localization.png`, `Figure_S_detection_curves_adult.png`, `figures_supplement.json` | Figs. S13–S15 |
| `report_figures_confirm.py` | The pre-specified spike run | `results/g4_confirm/g4_confirm_summary.json`, `g4c_<head>_summary.json`, `results/g4/g4_pediatric_comparison.json` | `results/report/Figure_R16_confirm.png`, `figures_confirm.json` | Figs. 10, S18 |
| `report_facts.py` | Merges the fact modules below into the named facts that `paper/build_paper.py` prints (value as printed, unrounded value, source) and checks that every source file exists; `--check`, `--json OUT` | the fact modules | none | every number in the paper |
| `report_facts_g12.py` | Facts of the spherical benchmark, the Hunold and Goldenholz adaptations, the adult comparison and regions, and the literature values | `results/g1a/`, `results/g1b/`, `results/g1c/`, `results/g2/`, `docs/literature/` | none | Sections 2–4; Table 3; S1–S3, S8, S9 |
| `report_facts_g3.py` | Facts of the head-size benchmark and the smaller heads | `results/g3a/`, `results/g3b/`, `configs/g2_adult.toml`, `docs/methods.md` | none | Sections 2–3; Tables 1, 4; S2, S4–S6, S10 |
| `report_facts_g4.py` | Facts of the exploratory spike study, localization and head motion | `results/g4/`, `docs/methods.md` | none | Section 3.6; Table 5; S7, S8 |
| `report_facts_methods.py` | Facts of the methods: noise model, OPM axes, time-domain model, resamples, seeds, parameter provenance | `results/`, `configs/`, `docs/methods.md`, `docs/provenance_register.md` | none | Section 2; S3, S4, S7, S8 |
| `report_facts_rev.py` | Facts of the noise-model checks and sensitivity analyses, the MRI check, the fitted helmet, the adult's depth-bin intervals and further seeds | `results/g2_noise_sensitivity/`, `results/g2_covariance_validation/`, `results/g3b_children_qc/`, `results/g3b_constant_gap/`, `results/g2/g2_depth_bins.json` | none | Sections 3.2–3.5; Tables 3, 4; S3–S6, S10 |
| `report_facts_confirm.py` | Facts of the pre-specified spike run | `results/g4_confirm/`, `configs/g4_confirmatory.toml` | none | Sections 2.5, 3.6; Table 5; S7.5, S10 |
| `report_facts_writer.py` | Further facts computed from stored values for the prose: the adult ratios in dB, counts in words, ranges over the principal heads and over the heads reported separately | `results/`, `configs/`, `docs/methods.md`, `docs/provenance_register.md` | none | throughout |
| `export_target_precision.py` | Full-precision depth below the scalp and bin membership of every row of the per-target tables | anatomies in `data/external/`, `cache/anatomy/`, the per-target tables | `results/g2/g2_targets_depth.csv`, `results/g3b/g3b_targets_<head>_depth.csv`, `results/g3b_constant_gap/targets_<head>_depth.csv` | Fig. 5 (through `study_g2_depth_bins.py`); depth counts |

## Publication

| Script | What it does | Main inputs | Outputs | Paper |
|---|---|---|---|---|
| `deploy_pages.py` | Builds the HTML version with `paper/build_html.py` from the rendered sources in `paper/build/` and commits it as the next commit of the branch `gh-pages` through a temporary index (the working tree, index and current branch are not touched); pushes only with `--push` | `paper/build/` (after `paper/build_paper.py`), `paper/figures/` | a commit on `gh-pages` | HTML version of the paper |
| `check_live_site.py` | Crawls the published HTML version and checks every internal page, figure, download and anchor (`--external`: external links too) | the published URL | report on standard output; exit status 1 on a failure | HTML version of the paper |
