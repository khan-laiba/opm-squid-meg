# Goal: Adult-to-pediatric OPM versus SQUID MEG comparison

## Mission and required sequence

I am extending my mentor Mainak Jas's work on on-scalp versus off-scalp MEG. Build a reproducible Python/MNE-Python study that first establishes a general adult comparison, drawing on Jas (2026), Hunold (2016), and Goldenholz (2009), and then extends that validated framework to pediatric anatomy.

The scientific sequence is:

**Adult analytical benchmark → realistic adult sensor/SNR comparison → pediatric fixed-helmet versus head-adaptive comparison → epilepsy-focused detection and bounded localization evaluation.**

The adult work is a substantive first result, not a preliminary illustration to skip. Pediatric anatomy acquisition may proceed in parallel, but do not substitute pediatric analyses for adult validation.

Determine how source depth, orientation, cortical extent, sensor geometry, and noise govern relative performance. Then test whether adapting sensors to a smaller head changes that performance relative to one fixed adult helmet.

Epilepsy, particularly interictal epileptiform discharges (IEDs), is the motivating application across both populations. Keep signal amplitude, sensor-level SNR, event detection, and source-localization accuracy distinct. Establish whether advantages exist, where they occur, and when they disappear or reverse; do not optimize the study to demonstrate OPM superiority.

Compare OPM magnetometers with Neuromag's 102 SQUID magnetometers and 204 SQUID planar gradiometers separately. Include the combined 306-channel SQUID array as a system-level comparator. EEG and axial gradiometers are outside scope.

## G0 — Audit the project and establish provenance

Read the three local papers in full, including equations, figure captions, and available supplements. Read the supplied TRIUX specification image and the hardware/implementation references below. Inspect existing repository instructions and any mentor-provided problem statement.

Audit `simulate_realistic_head.py`, `simulate_extended_sources.py`, existing anatomy, geometry, configurations, tests, and cached results. Reuse validated components. Preserve the earlier dense conformal arrays and constant-noise-ratio results as a labeled idealized baseline; numerical field-sampling locations are not independent physical sensor channels.

Maintain `GOAL.md`, a milestone-based `PLAN.md`, `docs/methods.md`, and a parameter/provenance register. Separate source-reported methods, hardware specifications, and new assumptions. Record unresolved ambiguities instead of silently repairing or combining incompatible methods.

**Acceptance:** An executable work plan identifies reference reproductions, adaptations, new experiments, required inputs, and numerical checks. Previous outputs are preserved, and unverified legacy claims are labeled as such.

## G1 — Establish the adult foundations from all three papers

### A. Mainak Jas: adult depth, distance, and noise benchmark

Reproduce the adult analytical results underlying Figures 3–4 using the supplied paper's definitions:

- Head radius 95 mm; brain radius 80 mm; tangential dipole 30 nAm.
- Ideal radial point magnetometers at 0-mm and 18-mm scalp standoff.
- Relative RMS noise ratio eta from 1 to 6.
- Source depth measured from the scalp, not from the sensors.

Independently compare analytical field maxima with numerical forward calculations. Verify the approximately 28-mm equal-SNR depth at eta = 3 to the paper's reporting precision. Exclude silent sources, including the sphere-center solution, from meaningful ratio/crossover calculations.

Document any assumed absolute noise needed to reproduce SNR-axis values; eta alone does not specify it. Retain the toy target-versus-background depth experiment as a separate explanatory benchmark, resolving text/caption differences explicitly.

Do not treat this idealized point-magnetometer calculation as a Neuromag gradiometer model. Reproduce experimental SEF results only if the actual authorized recordings and necessary metadata are available; otherwise use the experiment as context, not simulated validation data.

### B. Hunold: adult depth–orientation–spike simulations

Implement the MEG portions of the adult study: cortical-normal focal sources, small orientation-constrained patches, spike-like waveforms, and distributed background brain activity.

Preserve its definitions: source depth relative to the scalp; orientation relative to the nearby inner-skull normal; approximately 80-ms spikes; approximately 20-mm² patches, not 20-mm-radius patches; and its depth/orientation binning where anatomy supplies enough sources.

Keep the 600-nAm focal-source reference and the paper's patch-strength/background conventions in a dedicated reproduction configuration. Do not interpret them as universal physiological amplitudes.

First verify the magnetometer-versus-planar-gradiometer comparison; then add OPMs as an explicit extension.

Implement Hunold's channel selection and Hilbert-envelope SNR faithfully: choose the channel with the largest noise-free spike amplitude, then evaluate its preceding baseline. Its brain-noise-only reference omits intrinsic sensor noise. Its empirical 2.5 visual-detection threshold is not a universal event-detector threshold.

### C. Goldenholz: adult cortical SNR mapping

Implement focal and extended cortical-source maps using the paper's source and noise definitions. Preserve the 10-nAm focal reference and 10/16-mm geodesic-radius patches with 50 pAm/mm² moment density in a separate configuration.

Implement Equation 1's channel-averaged power SNR, including the `1/N` factor. Keep modeled background noise separate from recorded spontaneous-activity noise. Reproduce the recorded-noise branch only with suitable data and metadata.

Record each paper's head-model conductivities as printed; do not silently harmonize them. In particular, investigate rather than automatically replace Goldenholz's printed skull conductivity of 0.06 S/m.

Where original anatomy, recordings, or implementation details are unavailable, provide a documented adult methodological adaptation. Matching a trend on fsaverage is not exact reproduction of the original participants.

**Acceptance:** The adult analytical benchmark passes independent numerical checks; adult depth–orientation maps and focal/extended cortical maps run; paper-specific SNR definitions remain distinct; each output states its reproduction/adaptation status.

## G2 — Build the realistic adult OPM–Neuromag comparison

### Physical arrays

Use actual 3-D Neuromag positions and orientations from a suitable FIF or verified geometry file. Preserve the rigid helmet and valid device-to-head/head-to-MRI transforms. A 2-D plotting layout is insufficient. Label representative rather than installation-specific geometry honestly.

MRN specifies T3 magnetometer type 3024 and planar-gradiometer type 3014. Verify imported metadata and use the pinned MNE version's finite coil integration points, weights, orientations, and baseline. Do not relabel magnetometers as gradiometers.

Use the supplied brochure's typical white-noise amplitude spectral densities as one hardware reference scenario:

    Magnetometers: 3.5 fT/sqrt(Hz) = 3.5e-15 T/sqrt(Hz)
    Gradiometers:  3.6 fT/(cm sqrt(Hz)) = 3.6e-13 T/(m sqrt(Hz))

Keep typical values, guaranteed bounds, and measured spectra distinct. The brochure's 18-mm pickup-coil–to–room-temperature-surface spacing is not a uniform scalp-to-sensor distance. Compute the actual gaps and document the Dewar-clearance model.

Construct an explicit finite OPM array. Start with a 102-site scalp-normal matched-site configuration where feasible. Use verified geometry or a labeled assumption based on Jas's experimental approximately 7-mm sensing-center distance from the helmet's inner surface and 10-mm cubic cell. State any additional scalp gap and test finite-volume averaging.

Verify package clearance/packing or label physical feasibility unresolved. Use device-specific noise specifications when verified; otherwise run a declared sensitivity sweep, such as 7–30 fT/sqrt(Hz), without presenting it as one device's measured specification.

Separate matched-site/coverage controls, feasible channel-budget controls, and full-system comparisons. Report sites, measured axes, channels, and retained rank. Keep zero-standoff point sensors only in the idealized benchmark.

### Common sources and explicit noise

Use identical target and background source realizations across arrays within each anatomy and condition:

    y_m(t) = L_target,m q_target(t)
             + L_bg,m q_bg(t)
             + E_m e(t)
             + epsilon_m(t)

Evaluate intrinsic sensor noise, brain activity plus intrinsic noise, and their combination with environmental interference. Retain paper-specific brain-noise-only references separately.

Project common background sources through each array's forward model. Do not force equal sensor RMS, equal sensor SNR, or eta = 3 in the realistic comparison. Include an independent cortical-background baseline and a bounded spatially correlated-background extension.

Define source-area/covariance scaling so mesh refinement does not create extra background power. Preserve signed summation and cancellation within cortical patches. Compare fixed-total scalar moment and fixed moment density as separate source-strength conventions.

Use common analysis filters and bandwidth for the primary comparison, with sensor response and frequency-band sensitivity analyses documented separately. Integrate one-sided noise PSD times the squared composite filter response to obtain variance, and verify the simulated variance.

Apply environmental fields through the actual sensor response. Any preprocessing must transform signal and noise consistently, with attenuation and rank changes reported.

### Adult endpoints

Generate signal-versus-depth curves, depth–orientation heatmaps, cortical SNR maps, and relative-performance maps for focal and extended sources. Include actual adult helmet fit and a bounded, source-blind head-position sensitivity analysis.

Keep raw field amplitudes in T and gradients in T/m separate. In addition to the paper-specific metrics, implement clearly labeled general measures:

    Peak-channel RMS SNR = max_i |s_i| / sqrt(C_ii)
    Mean-power SNR_dB = 10 log10[(1/N) sum_i s_i²/C_ii]
    Known-topography detectability = sqrt(sᵀ C⁻¹ s)

Here s is a declared noise-free target topography and C is the matching noise covariance. Use rank-aware whitening, not unstable explicit inversion. Distinguish oracle and independently estimated covariance. The final expression is not calibrated event sensitivity or localization accuracy.

For combined SQUID data, retain cross-type covariance and appropriate numerical scaling. Verify dimensionless results are invariant to consistent unit changes; never average unscaled T and T/m.

**Acceptance:** A substantive adult report exists for every sensor configuration, with validated geometry/units/noise, numerical convergence checks, uncertainty, and explanations of differences from the idealized benchmarks. Freeze and version this baseline before pediatric outcome comparisons.

## G3 — Extend the same framework to children

### Separate head-size effects from fixed-helmet mismatch

First reproduce the head-size extension already present in Jas's Table 1/Figure 5, after validating the adult case. Its SQUID shell follows `s = h + 18 mm` as head radius h changes. Label this the size-following, constant-standoff benchmark.

Then run the new hardware experiment: keep the same adult Neuromag helmet, sensor dimensions, and intrinsic-noise assumptions fixed across anatomies, while fitting the OPM array to each head.

Do not shrink or conform the real SQUID helmet to a child. Conversely, do not shrink OPM cell/package dimensions to make a dense array fit.

Use age-appropriate, physically scaled anatomy, beginning with a school-aged model and adding younger/larger heads as validated data allow. Preserve native dimensions when computing forward models. An adult surface scaled down is a size-only control, not validated pediatric anatomy; a standard-space brain template is not automatically a full-head fit model.

Record actual head dimensions and anatomical limitations, not age labels alone. Obtain more than one pediatric anatomy for any claim about anatomical variability; template-only findings remain conditional simulations rather than population estimates.

### Model regional fit and positioning fairly

Define feasible neutral and bounded translated/rotated SQUID placements using physical clearance and available evidence. Include a reasonably well-fitted position, not only unfavorable placements. Select the positioning rule without knowledge of each test source.

Compute regional scalp-to-sensor gaps, source-to-sensor distances, coverage, and sensitivity. Keep these separate from source depth and orientation. Report which regions gain or lose sensitivity when the head moves nearer one part of the helmet.

Refit OPMs using a consistent coverage policy. Document any changes in achievable channel count, packing, or scalp coverage. Test their effects separately rather than attributing all differences to head size.

Use geometry-only controls with common source/noise conventions first. Add age/state-dependent background physiology only as a separate, evidence-supported sensitivity analysis. Do not alter several factors simultaneously and call the entire change a helmet-fit effect.

Where useful, add a labeled counterfactual array with reduced head–helmet mismatch while holding sensor-response/noise conventions constant. It is a mechanistic control, not a manufactured pediatric SQUID system. Distinguish this controlled contrast from the combined technology/geometry differences in real arrays.

### Test both pediatric performance and additional benefit over adults

Reuse the adult source families, analysis definitions, calibration policies, and reporting code. Match comparisons across anatomies using homologous regions and declared depth/orientation/patch-size strata, with native-space calculations and area weighting. Show unmatched or sparse strata rather than inventing correspondence.

For the same chosen dB SNR metric, evaluate separately for each SQUID comparator:

    D_child = SNR_OPM,child - SNR_SQUID,child
    D_adult = SNR_OPM,adult - SNR_SQUID,adult
    Delta   = D_child - D_adult

Report D_child, D_adult, and Delta, not Delta alone. A positive Delta indicates an increase in relative OPM performance under the specified matching assumptions; it does not necessarily mean OPM outperforms SQUID in the child.

Map where OPM, SQUID, or neither achieves useful absolute performance. Report uncertainty, noise/pose dependence, and model-conditional interpretation. Do not use ratios of dB values or turn a depth percentage into a cortical-area/brain-volume fraction.

Equal-performance boundaries may be regional, absent, or multiple. Do not impose a universal crossover depth or confuse a depth-binned summary with an individual-source crossing.

**Acceptance:** The validated adult framework runs on pediatric anatomy, and results distinguish head size, fixed-helmet fit, coverage, and background assumptions. Advantage, no detectable advantage, and disadvantage are equally valid outcomes.

## G4 — Evaluate epilepsy relevance without replacing the core study

Use the adult and pediatric forward/noise models to assess individual IED-like events across representative regions, depths, orientations, and patch extents. Sweep source strength and morphology rather than choosing only easily visible spikes.

Keep a known-source/onset oracle distinct from a practical detector that does not receive the true event time or source. Calibrate thresholds on independent null data and freeze them before held-out evaluation. Account for searches across time, channels, and templates.

Report sensitivity versus false events per minute using the same operating-point policy across systems and ages. Distinguish event-level false alarms from sample-wise false-positive probabilities. Include enough simulated duration/events to support uncertainty estimates.

An optional minimum-detectable-source-strength analysis may estimate the strength required for a prespecified detection probability, but any numerical threshold is an operational study choice, not an established clinical standard. Do not let a cortex-wide detector project delay the core SNR comparison.

On a bounded subset, evaluate focal equivalent-current-dipole localization and distributed MNE/dSPM reconstruction using the MNE examples as implementation references. Include off-grid truth and bounded registration/model mismatch; do not equalize sensor SNR before comparing accuracy.

Report localization error, patch-support recovery, failed fits, and nondetections. Distinguish localization among detected events from joint detection-and-localization success. Simulated IED-source recovery does not identify an epileptogenic zone or establish surgical benefit.

Study static fit first. A head-mounted array's motion advantage requires a separate time-varying analysis: OPM sensor–head geometry may remain fixed while its relationship to residual room fields changes. Motion robustness and sensor slippage are bounded secondary extensions, not conclusions established by static proximity maps.

**Acceptance:** Adult and pediatric epilepsy examples use the same validated framework; detection and reconstruction claims have separate supporting results and limitations.

## G5 — Reproduce, report, and prepare publication

Continue in the intended repository under `khan-laiba`, verifying the remote and private visibility. Create a new private project only if no appropriate one exists and access is verified. Preserve unrelated work and uncommitted changes; do not force-push or change visibility automatically.

Deliver modular scripts, pinned dependencies, configuration files, tests, provenance, numerical results, figures, and exact setup/smoke/full-run commands. Cache forward solutions using anatomy, geometry, transforms, coils, sources, and version hashes. Keep expensive runs resumable and separate from routine CI.

Test analytical fields, channel counts/types, transforms, physical clearance, uniform-field gradiometer cancellation, known-gradient response, source/background normalization, PSD/covariance, unit invariance, and numerical convergence. Use fixed seeds and repeated realizations. Repeated events and cortical vertices are not independent participants.

Build the project GitHub Pages report from validated result artifacts, in this order: adult reference benchmarks; realistic adult results; pediatric fixed-versus-adaptive extension; epilepsy implications; methods, uncertainty, limitations, and reproduction instructions.

Include downloadable numerical results and the generating commit/run identifiers. Verify site navigation, assets, captions, and downloads. Keep restricted recordings/anatomy, reference PDFs, secrets, and identifying metadata out of Git/public pages; verify redistribution permissions for derived assets.

Keep development private. Build locally or as access-controlled artifacts until release is authorized. Verify Pages visibility rather than assuming a private repository creates a private website. Public repository conversion and Pages deployment require explicit owner approval; record release-ready and publicly deployed as different statuses.

Finish with milestone evidence, repository/commit, exact commands, results/site locations, conditional scientific findings, and any failed, not-run, or blocked items. Do not mark a scaffold, attractive plots, or a pediatric-only analysis as completion.

**The first scientific deliverable is the general adult comparison. The final contribution is its controlled extension to pediatric fixed-helmet versus head-adaptive recording, with epilepsy relevance evaluated rather than presumed.**

## Source register

Source-reported methods and new study choices must remain distinguishable in the implementation and report.

- Jas et al. (2026), supplied August 21 preprint: *Signal-to-noise ratio of event-related fields in on-scalp and off-scalp MEG*. DOI: `10.64898/2026.08.17.744953`. Adult/head-size geometry and analytical SNR: PDF pp. 6–9, Table 1, Figures 3–5; experimental OPM geometry: p. 10; gradiometer and brain-noise extension discussion: pp. 23–25.
- Hunold et al. (2016): *EEG and MEG: sensitivity to epileptic spike activity as function of source orientation and depth*. DOI: `10.1088/0967-3334/37/7/1146`. Adult anatomy, sources, and noise: published pp. 1148–1151; SNR: section 2.3; limitations and comparison of SNR definitions: pp. 1158–1160.
- Goldenholz et al. (2009; online 2008): *Mapping the signal-to-noise-ratios of cortical sources in magnetoencephalography and electroencephalography*. DOI: `10.1002/hbm.20571`. Adult source/head model: p. 1078; SNR and noise: pp. 1079–1080, Equations 1–3; metric/positioning limitations: pp. 1081–1083.
- User-supplied Elekta Neuromag TRIUX technical-specifications image: typical white-noise values, distances, and array specifications. Record model/version applicability separately from other Neuromag generations.
- MRN Neuromag hardware: https://www.mrn.org/collaborate/elekta-neuromag-meg
- MNE implementation, units, and coil definitions: https://mne.tools/stable/documentation/implementation.html
- MNE inverse examples: https://mne.tools/stable/auto_examples/inverse/index.html
- GitHub Pages visibility: https://docs.github.com/en/pages/getting-started-with-github-pages/changing-the-visibility-of-your-github-pages-site

Before implementation, record the versions of external specifications/documentation actually used. This prompt defines a proposed study and does not claim that Mainak has reviewed or approved it.
