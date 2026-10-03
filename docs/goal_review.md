# End-to-end goal review (2026-10-02)

Review of the study against every milestone of `GOAL.md` (G0-G5) at commit 9771f68 of
`khan-laiba/opm-squid-meg` (private; tag `adult-baseline-v3`), written as the closing statement the
goal asks for in G5: milestone evidence, repository and commit, exact commands, result and report
locations, conditional scientific findings, and the items that failed, were not run or are blocked.
Nothing in this document claims that any author of the reproduced papers has reviewed or approved
the study.

How the review was done: (1) every requirement sentence of `GOAL.md` was checked against the code,
configuration, documentation and result files; (2) three independent reviewer agents audited G0-G1,
G2-G3 and G4-G5 with the cross-cutting rules, each in a fresh clone, read-only; (3) the test suite
was run in a clean clone of the commit; (4) the two complete external reviews of 2026-10-01 (Fable 5.1:
approve with notes; GPT 6 Astra via Codex: reject pending revision, then a re-check of the revision)
and the two independent verification passes of the revised documentation are taken into account
(`PLAN.md`, "Final reviews"). The gaps found here were fixed in the commit that adds this document
where the fix is documentation or a command; the rest are listed as open.

## Verdict per milestone

| Milestone | Verdict | One-line basis |
|---|---|---|
| G0 audit, provenance, plan | met | `docs/audit.md`, `PLAN.md`, `docs/provenance_register.md`, `docs/methods.md`; legacy outputs preserved and labelled idealized; ambiguities recorded (U-* rows) |
| G1A Jas analytical benchmark | met (REPRO) | Eq. 1 vs Sarvas 4.7e-15, vs MNE sphere 5.4e-8; d_eq(eta = 3) = 27.665 mm; sphere centre excluded; sigma recovered and documented (R-sigma); Fig. 6 toy benchmark kept with the caption difference recorded (U-J6); no SEF recordings, used as context only |
| G1B Hunold | met (ADAPT; OPM NEW) | sample anatomy, 20-mm^2 patches, 600-nAm reference in `configs/hunold_reference.toml`, Hilbert-envelope SNR with the max-noise-free-channel rule, MM vs GM first, OPM as an extension; calibration ambiguity recorded (U-HU-*) |
| G1C Goldenholz | met (ADAPT; OPM NEW) | Eq. 1 with 1/N, 10-nAm focal and 10/16-mm patches at 50 pAm/mm^2 in `configs/goldenholz_reference.toml`, 0.006 and 0.06 S/m both run and labelled, modelled vs recorded noise separate |
| G2 realistic adult comparison | met with reservations (NEW) | real Neuromag geometry and transforms, T3 coil types verified, brochure noise as one scenario, explicit finite OPM arrays with clearance verified, three noise conditions, three general metrics with rank-aware whitening, convergence and uncertainty, a report per configuration now covering all three metrics; frozen as `adult-baseline-v1/v2/v3`; reservations: the plug-in covariances are not paired across arrays, the 204-site control is practically the full system, no measured-spectrum SQUID scenario |
| G3A size benchmark | met (REPRO + NEW contrast) | Jas Table 1 / Fig. 5 reproduced; fixed-shell contrast labelled |
| G3B fixed helmet vs head-adaptive OPM | partially met (NEW) | adult helmet and noise fixed, OPMs refitted, native dimensions, size-only controls labelled, three templates of one database (no population claim), source-blind placements incl. a well-fitted one and a counterfactual helmet, D_child/D_adult/Delta per comparator with strata, usefulness maps; two reservations: the clearance rule left the adult's OPM sensors about 0.8 mm farther from the scalp than the children's (bounded by a sensitivity, below), and no native school-aged anatomy (owner decision, open item 1) |
| G4 epilepsy (adult, pediatric, motion) | met (NEW) | oracle vs practical detector, frozen thresholds, ROC per minute, S50 as an operational choice, bounded ECD/dSPM with mismatch, nondetections and joint success reported, no fit failed, motion as a bounded secondary extension; same framework for both ages |
| G5 software, reproduction, report | met for a private deliverable; release blocked on owner decisions | pinned environment, tests (133, clean clone 132 + 1 skip), caches keyed by content hashes, resumable runs, local report in the prescribed order with downloads and commits, private and no Pages (verified via the API and the documentation); no routine CI workflow (open) |

The scientific sequence the goal requires was followed: the adult comparison was established and
frozen (`adult-baseline-v1`, then v2) before any pediatric outcome was compared; v3 recomputed the
adult and the pediatric results from one code commit after the final reviews found faults in the
arrays and the detector scoring (not after any pediatric outcome); the pediatric analyses reuse the
adult framework unchanged; epilepsy relevance was evaluated, not presumed. The study was not
optimised for OPM superiority: the conditions under which the OPM does not win (intrinsic noise
only, 30 fT/sqrt(Hz), scalp gaps, the projected condition, the matched-site array at depth, the
counterfactual helmet, deep spikes, the adult localization) are reported next to those where it does.

## Milestone evidence

All result JSON files record the code commit that produced them (`provenance.commit`; a
`+dirty` suffix marks uncommitted code, present only in the two historic v1 diagnoses). The
OPM-dependent results come from commit 19a8fd2 (the G4 summaries re-summarised from the stored
simulations at d667a9d); G1A (5b48602), G3A (663eaa1) and the exact-sphere BEM check (81168f3) involve
no array and are unchanged.

| Milestone | Code | Results (`results/`) | Report pages (`site/_build/`) |
|---|---|---|---|
| G0 | `docs/audit.md`, `PLAN.md`, `docs/provenance_register.md`, `legacy/` (preserved) | - | `register.html` |
| G1A | `src/opmsquid/sphere.py`, `scripts/g1a_jas_benchmark.py`, `legacy/replicate_figure3.py` | `g1a/` (7 files: 4 figures, benchmark JSON, curves CSV, Fig. 3 replica folder) | `benchmarks.html` |
| G1B | `src/opmsquid/hunold.py`, `scripts/g1b_hunold.py`, `configs/hunold_reference.toml` | `g1b/` (19 files: 17 figures, summary JSON, sources CSV) | `benchmarks.html` |
| G1C | `src/opmsquid/goldenholz.py`, `scripts/g1c_goldenholz.py`, `configs/goldenholz_reference.toml` | `g1c/` (5 files) | `benchmarks.html` |
| G2 | `src/opmsquid/{g2,opm,neuromag,noise,background,environment,metrics,noisemodel,forward,fullres}.py`, `scripts/g2_*.py`, `scripts/study_*.py`, `configs/g2_adult.toml` | `g2/` (24 files: 14 figures, 7 JSON incl. BEM studies, 2 CSV, `G2_report.md`) | `adult.html`, `g2-report.html` |
| G3A | `scripts/g3a_jas_size_benchmark.py` | `g3a/` (4 files) | `pediatric.html` |
| G3B | `src/opmsquid/pediatric.py`, `scripts/g3b_pediatric_helmet.py`, `configs/g3b_pediatric.toml` | `g3b/` (21 files: 13 figures, summary JSON, 6 target CSVs, `G3B_report.md`) | `pediatric.html`, `g3b-report.html` |
| G4 | `src/opmsquid/{ied,detection,localization,motion}.py`, `scripts/g4_*.py`, `scripts/study_g4_*.py`, `configs/g4_epilepsy.toml`, `configs/g4_motion.toml` | `g4/` (51 files: 19 figures, 16 JSON, 13 CSV, 3 reports) | `epilepsy.html`, `motion-report.html` |
| G5 | `tests/` (133 tests), `scripts/run_all.sh`, `scripts/build_site.py`, `site/`, `requirements.txt`, `docs/release_checklist.md` | - | `index.html`, `methods.html`, `reproduce.html` (60 downloads with size, SHA-256 and commit) |

Exact commands (`README.md`, "Reproduce"): setup (`python3.12 -m venv .venv`, `pip install -r
requirements.txt`, the MNE sample data and the three infant templates), the smoke run (unit tests,
then the G1A benchmark, about 15 min) and the full run (`bash scripts/run_all.sh`, about 8 h). The
local report: `.venv/bin/python scripts/build_site.py` (into `site/_build/`, git-ignored, link-checked
by `tests/test_site.py`). Clean-clone test run of 9771f68 (2026-10-02, fresh clone, the
development venv, an empty cache): 133 tests, 132 pass, 1 skipped (the stored-lead-field
fingerprint check, which needs the local cache); the report built and link-checked inside the suite.

## Conditional scientific findings

Conditional on one adult head (the MNE sample subject at its measured head position), an assumed
OPM noise of 15 fT/sqrt(Hz) (7-30 swept), no helmet-to-scalp gap beyond the 7-mm standoff (0-6 mm
swept), a 3-layer BEM, modelled brain noise calibrated on the gradiometers, average infant templates
of one database and scaled copies of the adult; details, intervals and limitations in `README.md`
and `docs/methods.md`:

- Adult (G2): the dense 205-site OPM array has 1.11x [1.08-1.14] the known-topography detectability
  of Neuromag's 306 channels (1.6x at 10-15 mm depth, 0.98-1.04x below 35 mm); after the 8-term
  external-field projection 1.05x [0.99-1.10], a tie within its interval; the matched 95-site array
  1.00x (0.89x projected). The advantage disappears at 30 fT/sqrt(Hz) (1.00x) and reverses with a
  6-mm gap (0.91x) or without brain noise (0.74x).
- Pediatric fit (G3B): with the child in the fixed adult helmet, Delta = +0.59 to +1.25 dB (five
  heads, D_adult +0.85 dB); in a helmet scaled with the head, centred, -0.21 to +0.01 dB: the
  relative gain is the fixed helmet's fit, not head size as such. The clearance rule left the
  adult's OPM sensors 0.78 mm farther from the scalp than the children's; with that equalised
  (sensitivity) the Deltas are +0.47 to +1.12 dB in the fixed helmet and -0.34 to -0.06 dB in the
  scaled helmet, the same reading. Both systems' absolute detectability rises in the smaller heads.
- Spikes (G4): at 10-20 mm the dense array needs about a quarter less strength for 50 % detection
  (ratio 1.29 [1.03-1.51], 8/2 locations, p = 0.04 uncorrected; 14/2 with the oracle), the same in
  every smaller head (1.33-1.61); no deeper difference is established; the matched array is behind at
  45-70 mm (1/9). Dipole errors are similar across arrays (coregistration-limited); a dSPM gain seen in
  earlier runs did not reproduce in the adult; across the six anatomies 29 of the 30 localization
  comparisons with p < 0.05 favour an OPM array (uncorrected; half for weak, mostly undetected sources).
- Motion (bounded): a 10-mm head displacement in the fixed helmet costs Neuromag 1.4-1.7 dB
  uncompensated; a 3-deg OPM cap slip 0.1-0.5 dB; un-modelled in-band rotation costs the OPM 1 dB at
  about 0.02 deg RMS in a 1-nT residual field without correction, 0.4 deg after a field projection
  with 1-deg/1-% calibration errors.

These are simulations with declared assumptions, not device measurements; G3B is a geometry and
helmet-fit experiment and G4 a bounded recovery of simulated sources, with no diagnostic,
epileptogenic-zone or surgical claim.

## Failed, not run or blocked

Failed: nothing. Every planned run completed; the one run that stopped (the adult localization
started before the lead fields were rebuilt, 2026-10-01) wrote nothing and was rerun.

Not run or not available:
1. No native school-aged anatomy (owner decision 2026-09-30): the school-age case is a scaled adult,
   labelled a size-only control. Any claim about anatomical variability would need individual children.
2. The TRIUX specification image was not on disk; its values were transcribed from `GOAL.md` and
   flagged (U-HW1).
3. No verified OPM device noise specification: a declared 7-30 fT/sqrt(Hz) sweep instead.
4. The Jas SEF recordings, the Hunold and Goldenholz participants and the recorded-noise branch of
   Goldenholz are unavailable: adaptations, labelled ADAPT.
5. No routine CI workflow: the tests run locally (`run_all.sh`, the smoke command); the expensive
   runs are resumable and separate. Adding a GitHub Actions workflow is an owner decision (minutes
   on a private repository).
6. A clean-environment smoke test with a fresh venv was done on 2026-09-30 (87 tests then); later
   verification used clean clones with the development venv.

Blocked (owner decisions, `docs/release_checklist.md`): licence; identifying author metadata in the
git history and earlier PDFs with licensed fonts (a public release needs a fresh export or a
rewritten history; no force-push or history rewrite has been done); redistribution of the
infant-template and fsaverage derivatives, the CC-BY figure raster, the digitised Hunold waveform,
the literature extractions and the goal text; two documents naming local folders; and, after all of
these, visibility and Pages. Release-ready: no. Publicly deployed: no.

## Gaps found in this review

Fixed in the commit that adds this document:
- An explicit smoke command in `README.md` (tests, then G1A).
- Failed dipole fits stated: none (`docs/methods.md`, section 9).
- The versions of the external documentation used, recorded in the register (DOC-MNE, DOC-PAGES)
  with the consequence for Pages: on a user account any Pages site would be public.

Fixed from the G4/G5 audit (below): failed fits counted under declared criteria; the pediatric
report's held-out range computed from the data and its censored intervals printed with their
finite end; the detector's truth-model candidates, the true-peak localization and the dSPM
implementation's deviation from MNE disclosed; threshold uncertainty, operational wording and
coverage limits stated; the report site's epilepsy page and status texts.

Fixed from the G0-G1 audit: the G1B and G1C drivers read their constants from the reference
configurations (G1B rerun: identical to the committed results to 7e-14 relative, a round-off
difference of `600 * 1e-9` against the literal `600e-9`; G1C rerun: identical); the two G1C
side-effects (patch truncation near the skull, the medial wall) measured in the run and
documented; stale configuration comments; the inferred Hunold numerator, the task-baseline nature
of the recorded noise, the Jas conductivity entry and the medial-wall dipoles recorded; the audit
line on repository instructions.

Fixed from the G2-G3 audit: the standoff sensitivity (`scripts/study_g3b_standoff.py`,
`results/g3b/G3B_standoff_report.md`) bounding the adult-only standoff; the metric dependence
(peak-channel SNR, where Neuromag is ahead; mean-power SNR) and the depth profile after the
projection in the G2 report, the G3B report, the README and the methods; the template Delta
estimator, the plug-in covariance draws, the 204-site control and the package footprint stated;
the G2 report's representative-geometry line and benchmark link; stale constants in code comments.

Open (minor):
- `results/g1a/fig3/figure3_summary.json`, written by the preserved legacy replica, carries no
  status label or commit of its own (the download list inherits G1A's commit, 5b48602).
- No CI workflow (above); `run_all.sh` has no mid-run resume (the lead-field chunks and the G4
  states are cached).
- The dSPM estimator is the study's own implementation; a rerun of the bounded localization with
  `mne.minimum_norm` would show whether the deviation moves any error or p-value.
- The localization study's practical-detector thresholds (10 min of null data) have no held-out
  check; localization against Neuromag combined only.

## Requirement audit

### G4 and G5 with the cross-cutting rules (independent reviewer 3; verdict: accept with reservations, G4 and G5 acceptance met)

| Requirement | Status | Evidence and disposition |
|---|---|---|
| IED-like events across regions, depths, orientations, extents; strength and morphology swept | met / partial | 1,728 events per anatomy, 10-320 nAm x 3 morphologies, focal and 10-mm patches; locations frontal-heavy, one patch extent: stated as a coverage limit (methods 9) |
| Oracle distinct from a practical detector without time or source | met | `detection.py`; the candidates' topographies come from the truth forward model: now disclosed (methods 9, A-G4-DET) |
| Thresholds on independent null data, frozen; searches over time, channels, templates accounted for | met | 10 + 20 + 20 min of fresh null draws; the per-sample maximum over templates x candidates; Poisson uncertainty of the thresholds now stated |
| Sensitivity vs false events per minute, same policy across systems and ages; event-level vs per-trial false alarms | met | ROC per anatomy (the five pediatric ROC figures now on the report site), held-out rates 0.40-1.55 per minute, matched-rate check |
| Minimum detectable strength as an operational choice | met | S50 and the 1-per-minute point labelled operational, not clinical (methods 9, A-G4-DET, the site) |
| Bounded ECD and dSPM with off-grid truth and mismatch, no SNR equalisation | met with a deviation | `mne.fit_dipole`; the dSPM estimator is the study's own implementation (depth weighting without MNE's limit): now documented as a deviation (methods 9, A-G4-LOC); the inverse at the true peak sample: now disclosed |
| Localization error, patch-support recovery, failed fits, nondetections | met (was partial) | failed fits were never counted; now counted under declared criteria from the stored events (`scripts/study_g4_fit_failures.py`, `results/g4/G4_fit_failures_report.md`): 24 of 576 dipoles per array more than 30 mm off among detected events |
| Localization among detected vs joint success; no epileptogenic-zone or surgical claim | met | separate medians and joint shares; the disclaimer now also on the epilepsy page |
| Static fit first; motion as a bounded secondary extension | met | methods 12, exact rigid-motion check |
| G5: repository, privacy, no force-push, no visibility change | met | `gh api`: private, no Pages; reflog without forced pushes |
| G5: modular scripts, pinned dependencies, configs, tests, provenance, results, figures; setup/smoke/full-run commands | met (smoke command added) | `requirements.txt` pinned; 133 tests; `README.md` |
| G5: caches keyed by hashes; resumable runs separate from CI | met / partial | forward cache and lead-field fingerprints keyed by geometry, transforms, coils, sources, versions; no mid-run resume of `run_all.sh`; no CI workflow (open) |
| G5: required test categories; fixed seeds; repeated realizations; events and vertices not independent participants | met | all ten categories present; the location and the parcel are the statistical units |
| G5: report in the prescribed order with downloads and commits, navigation verified | met | 11 pages, 62 downloads with SHA-256 and commit, 53 figures with captions, link check in the suite |
| G5: restricted data, PDFs, secrets, identifying metadata kept out; redistribution verified | partial (owner-blocked) | nothing restricted tracked; commit author metadata and derived-asset permissions are the documented release blockers |
| G5: private development, Pages verified, two statuses, the finish statement | met | `docs/release_checklist.md`; this document |
| Cross-cutting: sequence, distinct quantities, advantages and reversals, comparators, scope, labels, spec versions, no approval claim | met / partial | v1 and v2 frozen before the pediatric outcomes, v3 recomputed adult and pediatric from one commit after review findings about the arrays; G4 localization against Neuromag combined only (stated); TRIUX values applied to the Vectorview geometry with MRN T3 coils (now stated, HW-mag-noise); documentation versions now recorded (DOC-MNE, DOC-PAGES) |

Reviewer gaps and their disposition: (1, major) failed fits: fixed above; (2) the pediatric report's
hard-coded held-out range: now computed from the data; (3) dSPM not MNE's: documented as a
deviation (a rerun with `mne.minimum_norm` is listed as open); (4) the truth-model candidates and
the true-peak localization: disclosed; (5) threshold uncertainty: Poisson intervals stated; (6)
operational wording: added; (7) coverage: stated as limits; (8) the report site: pediatric ROC
figures, operating-point notes, the disclaimer, the status texts; (9) smoke command and the
failed/not-run/blocked list: this document and `README.md`; (10) documentation versions: recorded;
(11) notes: the censored interval format in the pediatric report fixed; no CI and the localization
detector's unchecked thresholds listed as open.

### G0 and G1 (independent reviewer 1; verdict: accept with reservations, G0 and G1 acceptance met)

The reviewer reproduced the G1A, G1B and G1C numbers independently (equal-SNR depths, sigma, the
Fig. 3 marker, the legacy cap d_eq values, the per-bin dipole counts, the depth and orientation
descriptors, the calibration ratio, the Eq. 1 median and share).

| Requirement | Status | Evidence and disposition |
|---|---|---|
| G0: papers read in full; legacy audited and reused; conformal results labelled idealized; GOAL/PLAN/methods/register; sources vs assumptions separated; ambiguities recorded; outputs preserved; legacy claims labelled | met | `docs/literature/*`, `docs/audit.md`, `legacy/README.md`, the register's [P]/[C]/[A] tags and U-rows |
| G0: the TRIUX image read | not met (blocked) | not on disk; values transcribed from the goal text and flagged (U-HW1); now also stated that they are applied to the Vectorview geometry with MRN T3 coils |
| G0: repository instructions and mentor statement inspected; documentation versions recorded | partial (now recorded) | a line in `docs/audit.md`; DOC-MNE and DOC-PAGES in the register; "before implementation" cannot be shown |
| G1A: parameters, standoffs, eta, depth from the scalp, analytical vs numerical, 27.665 mm, silent sources excluded, sigma documented, toy benchmark, no gradiometer claim, SEF as context, status on outputs | met | `opmsquid.sphere`, `results/g1a` |
| G1B: sources, patches of ~20 mm^2, binning, 600-nAm reference in a dedicated configuration, MM vs GM first, OPM as extension, channel selection and Hilbert SNR, brain-noise-only reference, 2.5 not a detector threshold, status on outputs | met with notes | the configuration is now the source of truth (the driver read some values from code); the inferred numerator stated in PLAN; patch bins not matched to the paper and sparse bins now stated (methods 6); the 171 medial-wall dipoles recorded (U-HU-wall) |
| G1C: paper definitions, 10/16-mm patches at 50 pAm/mm^2 in a configuration, Eq. 1 with 1/N, modelled vs recorded noise, conductivities as printed with 0.06 S/m investigated, status on outputs | met with two measured side-effects | the usable-vertex rule truncates patches near the skull and G1C keeps the medial wall that G2 and G4 exclude: both now computed in the run (`patch_truncation`, `without_medial_wall` in `g1c_summary.json`: 22-31 % of the patches lose more than 5 % of their area, 7-9 % of patch SNRs change by more than 1 dB with a median of 0; without the wall the focal median is -21.66 instead of -22.13 dB) and stated (methods 7, A-BEM-DIST); the recorded-noise branch named task-baseline noise (GO-rec); 0.006 S/m stated as the headline value (PLAN) |

Reviewer gaps and their disposition: m1 patch truncation and m2 medial wall: measured and
documented (a G1C variant without the wall and with untruncated patches is listed as open); m3
configuration not the source of truth: fixed; m4 TRIUX image: owner-blocked; m5 documentation
versions: recorded; m6 audit line: added; m7 status labels: the CSVs carry no status line (their
JSON does; open), 0.006 S/m stated; m8 numerator wording: fixed, the segment-only Hilbert variant
(about 0.6 %) listed as open; m9 recorded-noise suitability: documented; m10 sparse patch bins:
stated (a minimum-count rule listed as open). Notes: legacy originals in history (4023b16), the g0
tag after the first G1 code, the Jas conductivity entry (added, J-cond).

### G2 and G3 (independent reviewer 2; verdict: revise; G2 met with reservations, G3 partially met)

The reviewer reproduced the G2 headline (1.11x [1.08-1.14]; projected 1.05x [0.99-1.10]; matched
1.00x), the ENBW and in-band noise, the G3A benchmark and all 15 primary G3B rows.

| Requirement | Status | Evidence and disposition |
|---|---|---|
| G2: real geometry and transforms; coil types verified; brochure noise as one scenario; actual gaps; matched 102-site configuration; standoff, cell and gap as labelled assumptions; finite-volume averaging tested; declared noise sweep; three noise conditions; no forced equality; independent and correlated backgrounds; mesh-invariant scaling; signed patches; both strength conventions; common filter and variance check; environmental field through the real response with attenuation and rank; the endpoints; the three general metrics; cross-type covariance; unit invariance | met | `results/g2`, `tests/`, methods 8 |
| G2: clearance verified or feasibility labelled unresolved | met (wording fixed) | exact whole-cell clearance; the package footprint now stated as unverified in the report and the register |
| G2: identical source realizations across arrays | partial | oracle covariances share one source covariance; the plug-in covariances are estimated from separate draws per array: now stated as a limitation (methods 8) |
| G2: matched / channel-budget / full-system separated | partial | with 205 dense sites the 204-site control is practically the full system: now stated; a 306-channel control would need multi-axis OPMs (not modelled) |
| G2: a substantive report per configuration with all three metrics; differences from the idealized benchmark | met (was partial) | the G2 report now gives, per OPM configuration, detectability, peak-channel SNR (Neuromag ahead, 0.88x) and mean-power SNR (1.06x), the depth profile after the projection (behind Neuromag from 30 mm) and the link to the analytical benchmark; README and methods state the metric dependence and the deep reversal |
| G2: baseline frozen before pediatric comparisons | met with a caveat | v1 and v2 frozen before the pediatric outcomes; v3 recomputed both after review findings about the arrays (stated) |
| G3: Jas size benchmark; adult helmet and noise fixed; no shrinking; native dimensions; scaled adult labelled; dimensions recorded; more than one anatomy, no population claim; placements incl. a well-fitted one, source-blind; regional gaps and distances; coverage and channel-count effects; geometry-only first; counterfactual labelled; homologous strata, area weighting, sparse strata shown; usefulness maps; uncertainty and pose; no dB ratios | met | `results/g3b`, methods 10 |
| G3: the same OPM rules on every head | partial (measured) | the clearance rule moved the adult's sites outward (7.8 vs 7.0 mm median height), not the children's: D_adult carries about 0.8 mm more standoff; bounded by the standoff sensitivity (`results/g3b/G3B_standoff_report.md`): with the children's arrays moved out by 0.78 mm, D_child falls by 0.11-0.18 dB, the fixed-helmet Delta stays positive (+0.47 to +1.12 dB instead of +0.59 to +1.25) and the counterfactual-helmet Delta is negative in every head (-0.34 to -0.06 dB instead of -0.21 to +0.01): the reading of G3B does not change |
| G3: a school-aged model first | not met (owner-blocked, labelled) | the school-age case is a scaled adult, labelled a size-only control throughout |
| G3: D_child, D_adult, Delta per comparator | met (estimator now stated) | for the templates Delta is the adult-area-weighted median of parcel differences, not the difference of the two medians: stated in README and methods |
| G3: advantages and disadvantages equally prominent | met (was partial) | peak-channel SNR (Neuromag ahead) now in the G3B report; the standoff asymmetry stated |
| G3: no depth-percentage-to-volume conversion | partial | G3A reports a homogeneous-sphere volume share (labelled U-J4); kept as a sphere quantity, not a cortical share |

Reviewer gaps and their disposition: (1, major) the adult-only standoff: measured, bounded and
stated; a geometry-consistent rerun (outer surface rebuilt from the MRI scalp) is listed as open;
(2, major) reporting balance: fixed in the reports, README and methods; (3, major) school-aged
anatomy: owner-blocked; (4) template Delta estimator: stated; (5) plug-in draws: stated as a
limitation; (6) the 204-site control: stated; (7) "densest feasible": qualified; (8) SQUID
measured-spectrum scenario: open; (9) README standoff wording: fixed; (10) G3A volume share:
labelled; (11) G3B patch families (no 20-mm or fixed-density patches): open; (12) placement
evidence uncited, no yaw variant: open; (13) representative geometry in the G2 report: added;
(14) the G2 report's cross-references: added; (15) template averaging and depth mix: stated.
Notes: stale constants in code comments fixed; the 2.5-mm matched-site statement clarified.

## Addendum (2026-10-03): equal OPM standoff (`adult-baseline-v4`) and the open items

After this review the owner asked for the G3B standoff asymmetry to be fixed and every pending item
done. What changed, and what it changes in the verdicts above:

Equal standoff. The infant templates' BEM head surfaces are built on their MRI head surfaces (each
vertex a vertex of it); the sample subject's stored outer skin lay a median 0.8 mm outside its MRI
scalp over the sensor region. The whole-cell clearance rule therefore pushed the adult's OPM
sensors outward (median height above the MRI scalp 7.78 / 7.76 mm vs 7.00 mm in every child), and
the stored outer skin's vertex normals, which are not its triangles' normals (they lie close to
radial), tilted the adult's sensitive axes by a median 8-9 deg against the local scalp (the
templates' by under 1 deg). From v4 every anatomy's head surface has its vertices on its MRI scalp
(A-BEM-CONFORM; the identity for the templates, inherited by the scaled controls): every OPM array
sits at a median 6.99-7.01 mm and its axes within about 1 deg (median) of the scalp plane. The rule
is defined by geometry alone and was fixed before any v4 result existed, but it changes the adult
baseline after the pediatric outcomes were known: a deviation from the freeze sequence, recorded in
the decisions log. Every result that depends on the head model or the OPM arrays was recomputed
(phase A at ed852b2: lead fields, adult and pediatric detection, BEM studies, motion; phase B at
e53bea8: G1B, G1C, G2, G3B, every localization; the cross-reading summaries and a decomposition of
the change at 236068d/9512f1b).

Effect on the findings (v3 -> v4). G2: the dense array 1.11x -> 1.14x [1.12-1.17] Neuromag combined;
after the external-field projection 1.05x [0.99-1.10] (a tie) -> 1.12x [1.08-1.15]; the matched array
1.00x -> 1.01x (a tie), 0.89x -> 0.95x projected. The shallow ratios did not change (1.6x at 10-15 mm);
the deep ones rose. A decomposition (`scripts/study_head_surface_effect.py`) attributes the
projected-condition change mainly to the axes (nearly radial axes had made the room-field patterns
alike those of deep sources), partly to the sites (standoff, and 12 lower occipital sites the
inflated outer skin had excluded: about 2 %), and none of it to the conductor surface; the 1-layer
BEM, which has no head surface, moved the same way. G3B: D_adult +0.85 -> +1.00 dB, while the
templates' D did not change and the scaled controls' rose by 0.04-0.05 dB, so Delta fell by
0.11-0.27 dB to +0.44 (school-age size), +1.07 (2-year size), +0.73, +0.88 and +0.96 dB (24, 18 and
12 months), close to what the v3 standoff sensitivity predicted (+0.47 to +1.12 dB); after the
external-field projection it fell more (+1.00 to +1.81 -> +0.42 to +1.09 dB), as the adult's
projected D rose (G2). The counterfactual helmet about the laterally centred head now reverses the
gain in every child (-0.16 to -0.38 dB; the templates' intervals include 0), so the reading (the
gain is the fixed helmet's fit) is unchanged. With the yaw variants added, top contact is the adult's second-lowest
of 12 source-blind placements (v3: third-lowest of 10); against each head's family median the differences of the medians are
0.14-0.20 dB smaller. G4: the adult's
superficial detection advantage 8/2 -> 11/1 locations (p = 0.004), and in v4 also at 20-30 and 45-70
mm (uncorrected; of the earlier runs only an early v2 run found one, at 20-30 mm, so reported, not
established); the matched
array's deficit at 45-70 mm (1/9) is not seen (3/6); localization: the adult's dSPM errors
for 320-nAm patches are 5.6 and 2.4 mm smaller with the dense and matched arrays (p = 0.054 and 0.067;
with MNE's own dSPM 1.0 and 2.3 mm; v3 0 mm), still not robust; across the six anatomies the direction
holds (52 of 240 comparisons with p < 0.05, 51 favouring an OPM array; eight survive a correction
within their anatomy and array, two across the children). The templates' detection outcomes are essentially
those of v3 (their arrays did not change); the scaled controls' superficial strength ratios are 1.60 and
1.65 (v3 1.40 and 1.61).

Open items closed: G2's three reservations (plug-in covariances now from one common realization;
a measured-spectrum Neuromag noise scenario; a triaxial channel-count control at Neuromag's channel
count); G3B's standoff reservation; G1C variant without the medial wall and with untruncated
patches; G1B minimum-count rule and segment-only Hilbert variant; bounded localization with
`mne.minimum_norm`, with Neuromag magnetometers or gradiometers alone, and with a held-out check of
the detector thresholds; G3B 20-mm and fixed-density patches, a yaw placement variant and citations
for the placement; status lines in every result CSV and on the G1A replica's summary;
`RESUME=1 scripts/run_all.sh`. An independent code review of the v4 changes (no critical finding;
three major, eight minor) was addressed before the cross-reading summaries were run.

Verdicts after v4: G2 met (the reservations above addressed); G3B still partially met, now only for
the missing native school-aged anatomy (owner decision; a download would need approval); G5 met for
a private deliverable (148 tests; a clean clone of db3f0dc passes 147, one skip without the lead-field cache; no CI
workflow, an owner decision); the rest unchanged.

## Addendum (2026-10-03, later): a native school-aged anatomy

The one remaining G3 gap was the missing school-aged native anatomy. With the owner's approval three
typically developing children of OpenNeuro ds005234 (Fadeev et al. 2024; 7.8, 8.3 and 8.7 years) were
added as children A-C. Their own white surfaces, aparc labels and dense MRI scalp are used. The dataset
has two problems, both measured before use: the snapshot's file tree shifts each subject's watershed
BEM and talairach.xfm into the preceding subject's folder (the BEMs matched by the surfaces' volume
information), and the watershed inner skull lies just below the scalp (a median 0.8-2.5 mm over the
upper head, against 9.7 mm in the adult; `results/g3b/school_anatomy_checks.json`). The skull is
therefore modelled (A-BEM-CHILD: inner skull at 8 mm below the scalp where shallower, no vertex within
2 mm of a white-surface vertex), which on the adult, whose skull is segmented, changes the G2 headline
by at most 0.003; the fiducials are the adult's transferred by a cortex fit (A-G3-FID), which
reproduces the templates' own to 1.6-14.4 mm and their head frame to 3.2-5.9 deg. Adding the children
left every number of the six earlier anatomies unchanged (97,048 values).

Result: at a similar head circumference (520 and 535 mm against 525 mm) children A and C gain less than
the scaled school-age control: Delta +0.30 [+0.11, +0.38] and +0.17 [-0.07, +0.34] dB against +0.44
[+0.34, +0.56] dB (the intervals overlap or touch); child B, with a smaller head (486 mm), gains +0.41
[+0.28, +0.54] dB. Reweighted to the adult's depth mix the children's differences are -0.05 to +0.05 dB.
Their cortex is nearly adult-sized (1,666-1,869 cm^2 against 1,878 cm^2; scaled control 1,470 cm^2),
so with the background fixed per unit area Neuromag's brain noise stays near the adult's; the size-only
control, which shrinks the cortex with the head, gains more than the two children of its head size.
In the counterfactual helmet about the laterally centred head children A and B lose (-0.76 and -0.59
dB) and child C, whose helmet barely shrinks, is level (+0.04 [-0.20, +0.20] dB), so the earlier
statement that it reverses the gain in every child holds for A and B but not for C.
Pediatric G4 (same detectors, seeds and inverse as the other anatomies): at 10-20 mm more locations
favour the dense OPM in each child (11/4, 12/3 and 14/2; strength ratio 1.27 [0.93-1.63], 1.18
[1.02-1.51] and 1.43 [1.06-1.62]; p = 0.022, 0.0085 and 0.0037, none surviving the correction over 24
detection comparisons); there is no deeper difference with the practical detector; both arrays' point
estimates of the strength for 50 % detection are higher than the adult's (Neuromag at 20-30 mm 95-132
against 86 nAm, dense OPM 82-122 against 72 nAm), with overlapping intervals and no difference claimed.
Localization: dipole errors similar (4.4-6.7 mm for detected 320-nAm focal spikes); dense dSPM errors
2.9-4.7 mm smaller (p = 0.071-0.25; MNE's dSPM 4.9-5.5 mm, p = 0.009-0.020); 13 of the children's 120
localization comparisons have p < 0.05, all favouring an OPM array, none surviving a correction
(`docs/methods.md` section 11).

Verdict after this addendum: G3 met with reservations. The school-aged anatomy is native in its cortex
and scalp, but its skull is modelled and its fiducials transferred; three children of one dataset show
between-child differences without estimating a population.
