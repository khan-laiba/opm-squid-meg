# Report progress

**Goal.** Publish the pediatric epilepsy research report as a GitHub Pages site built from the
analysis already in this repository, as one coherent story that starts from Jas et al. (2026) and
builds step by step to pediatric epilepsy, following the story plan and the study questions
(the brain regions covered by OPM epilepsy studies, OPM vs SQUID by region across head sizes, the
number of sensors). Every number and figure must trace to an analysis output in the repository.

**Done when** the site builds and deploys, the live URL serves the final report with every link
and figure resolving, and two independent referees, Claude Fable 5.1 (maximum effort, a fresh
subagent) and GPT-6 Astra (maximum reasoning, Codex CLI), each return "accept" or "minor
revisions" in the same round. Each referee sees only the files.

**Site:** https://khan-laiba.github.io/opm-squid-meg/ (served from the `gh-pages` branch; the final
version replaces the round-1 version). The repository is public (code MIT; report,
figures and result files CC BY 4.0).

## Status

| Step | Status | Notes |
|---|---|---|
| 1. Foundations: number provenance, new figures from stored outputs, verified literature table, site integration | done | 31,150 sourced facts; figures R1-R10; 15 studies (14 verified in full text); 197 tests pass |
| 2. Draft the report as one narrative | done | about 8,500 words of main text, 10 figures, 3 tables, 39 references |
| 3. Internal verification: numbers, claims against evidence, methods completeness, build, tests, links | done | 166 checker issues resolved; supplementary pages cleaned for publication; 197 tests pass |
| 4. Referee rounds (both referees, same round) | done | round 1: major revisions from both; round 2: minor revisions (Fable 5.1), major revisions (GPT-6 Astra); round 3: minor revisions from both; round 4, on the deployed final version (f167fd3): minor revisions from both |
| 5. Deploy to GitHub Pages and check the live site | done | final version: acb4752 (round-4 corrections applied) deployed as gh-pages f8292c3 and checked live: 14 pages, 241 internal assets, 139 anchors, 0 failures |

## Referee rounds

- Round 1 (commit f59f76e): **major revisions** from both referees. They verified the numbers (no
  number contradicted its result file) and asked for: evidence on the noise model that decides the
  comparison (validation against measured Neuromag noise; near-skull cortex, coloured OPM noise and
  cardiac/ocular fields as sensitivity analyses); a quality check of the school-aged children's MRI
  surfaces; a helmet fitted to each head at the adult's gap, with placement uncertainty; a confirmatory
  spike run with the endpoint fixed in advance; public code; and a shorter manuscript in its own terms.
- Round 2 (commit 1d7230c): **minor revisions** (Claude Fable 5.1) and **major revisions** (GPT-6 Astra).
  Both again found no number that contradicts its result file and judged the first-round points largely
  addressed. GPT-6 Astra's two major points: the school-aged children's anatomy cannot carry a mechanism
  (measured to the MRI head boundary, child B's shallow targets fall from 247 to 47), and the helmet
  fitted at the adult's gap shows that fitting reduces the smaller heads' extra advantage, not that the
  advantage needs the wide gap. Fable 5.1 asked for multi-axis clinical arrays to be named, a checkable
  pre-specification, the excluded cortex per head and a plainer Abstract.
- Round 3 (commit fabd710, after two internal pre-checks): **minor revisions** from both referees (Claude Fable 5.1
  and GPT-6 Astra), the acceptance condition. Both again found every checked number in its result file (Fable 5.1
  recomputed about 30 groups of values, GPT-6 Astra 29) and judged every second-round point adequately addressed. Their
  remaining points are presentational. The small wording corrections applied before deployment: the Ren et al.
  comparison no longer read as selective compatibility, the template-convention bias stated as possible rather than
  certain, the 4,680-sample covariance described as the declared 60-s scenario, a pointer from Figure 7 to the
  interaction in Table S6, child C's within-head contrast attributed to the construction, the children's depth-strata
  count, and 'an independent prediction' for the magnetometer brain noise (fb7fd2f). The referees' remaining minor
  points were then applied in a final pass (log, 2026-10-05), and the final version went to both referees once more
  as a confirmation round.
- Round 4 (commit f167fd3, the deployed final version): **minor revisions** from both referees (Claude Fable 5.1 and
  GPT-6 Astra), again the acceptance condition, with every third-round point judged answered. Their remaining requests
  are editorial: the Abstract within NeuroImage's 250 words, the 18-month template's MRI-check flag beside its results,
  the children's skull rule as implemented, Section 3.6's checks one at a time with each array's false-event rates, the
  exact-p file named as the primary source of Table 3's p values, the superseded clinical-comparison record marked, and
  items only the author can supply (funding, competing interests, contributions, an archive DOI, print figure files).

## Log

- 2026-10-04: Pages enabled on the repository (served publicly, repository private); Codex CLI
  installed and logged in; the publishable files checked for local paths and e-mail addresses (none);
  shared interfaces added (`scripts/report_facts.py`, `scripts/report_style.py`).
- 2026-10-04: step 1 started as eight parallel agents, each on its own files: three fact modules
  (`scripts/report_facts_g12.py`, `_g3.py`, `_g4.py`: every number with its result file and key),
  two figure scripts (`scripts/report_figures_adult.py`, `_pediatric.py`: new figures drawn from
  stored outputs into `results/report/`), the verified literature table
  (`docs/literature/epilepsy_opm_studies.md`), the report page in the site builder (the report as
  the landing page, the earlier pages as supplementary material, a number-provenance page), and
  the deployment and live-site check tools. Next: one writer drafts the report, nine independent
  checkers audit it (number provenance, claims, physics, methods, references, flow, publication
  safety), one editor revises; then the referee rounds.
- 2026-10-04: step 1 done. The fact modules recomputed every number of the private plan files from the result
  files and corrected a few of them (e.g. the adult's placement spread is 0.46 dB, not 0.45); ten new figures
  (results/report/Figure_R1-R10) each checked against the result files; the literature table covers 12 clinical
  OPM epilepsy studies and 3 modelling precedents. All 197 tests pass.
- 2026-10-04: steps 2 and 3 done. One writer drafted the report from the facts, figures and verified literature; nine
  independent checkers (number provenance x3, claims and statistics, MEG physics, methods, references, flow,
  publication safety) raised 166 issues, which an editor resolved (two critical contradictions about the pediatric
  conclusion and the spike result among them). The abstract was shortened, the Introduction now opens from
  Jas et al. (2026), and the supplementary pages got licence and data credits and lost internal process notes.
- 2026-10-04: first deployment. The site built from f59f76e was committed to `gh-pages` (33c192f) and Pages
  set to serve that branch. Live check: 13 pages, 157 internal assets and 95 anchors, 0 failures; 42 external
  links, 32 load from a script and 10 publisher pages (Wiley, AIP, RSNA, MDPI) refuse scripts (to be checked in a
  browser). Round 1 sent to both referees with the manuscript, its site and the analysis files only.
- 2026-10-04: round 1 returned: major revisions from both referees (reports kept outside the repository).
  The owner made the repository public; code under the MIT licence, report, figures and result files
  under CC BY 4.0 (LICENSE, LICENSE-CONTENT.md).
- 2026-10-04: the analyses the referees asked for were implemented with tests (10e37b9): the children's MRI
  quality check (their T1 and head-mask files fetched after the owner's approval, with SHA-256 in
  configs/school_subjects_qc_manifest.json), the noise model against the measured Neuromag covariance, the
  noise-model sensitivity analyses, a helmet fitted to each head at the adult's gap, and a confirmatory spike
  run whose endpoint, seeds and number of locations were declared in configs/g4_confirmatory.toml before it ran.
- 2026-10-05: results of the checks committed (ae458a9). The children's surfaces are usable and their shallow
  cortex is genuine; Neuromag's detectability with its measured noise is 1.14 times the model's, so the adult
  OPM advantage is a conditional prediction whose size matches the model's own error for Neuromag; the named
  omissions (near-skull cortex, coloured OPM noise, heart and eyes) barely move it. The helmet fitted at the
  adult's gap brings the smaller heads' OPM advantage back to the adult's (and below it in two of the three
  children). The confirmatory spike run is computing (nine anatomies, one clean commit).
- 2026-10-05: revision of the manuscript started on branch `revision-r2`: clean main-text figures, the
  Methods details as sourced facts, the clinical comparison as qualitative analogues, a supplementary text
  page, and the condensed manuscript (about half the length), followed by independent checks before round 2.
- 2026-10-05: the confirmatory spike run finished (nine anatomies at one clean commit, ae458a9, settings as
  declared in configs/g4_confirmatory.toml). The declared endpoint (dense OPM array vs Neuromag's 306 channels,
  practical detector at 1 false event per minute, spikes 10-20 mm deep) passes the Holm-corrected test in all
  nine anatomies, with S50 ratios of 1.27 to 1.47, and in each of five independent noise replicates; the
  oracle and a mismatched detector agree; the site-matched array passes in five of nine. Results committed
  (e9433e5) with those of the helmet fitted at the adult's gap.
- 2026-10-05: stage B2 of the manuscript revision done: the manuscript rewritten at about half its length
  (about 5,500 words) with the supplementary text S1; nine independent checkers (numbers x3, claims, referee
  coverage, figures and tables, flow and style, compliance, build) raised 13 blockers and 71 major issues,
  an editor resolved them, three rechecks and a second edit followed. New figures from stored outputs: the
  arrays (R15), the geometry with the fitted helmet (R12), the noise-model checks (R17), the confirmatory
  spikes (R16), the children's MRI check and three supplementary spike figures; dB axes on the adult figures.
  The 3-70 Hz band was added to the band analysis and the children's MRI check re-run with full-precision
  depths, both from a clean checkout. Stage B3 (the confirmatory results written in, nine checkers, editor,
  point-by-point response letter) is running; round 2 follows.
- 2026-10-05: stage B3 done: the confirmatory results written into Section 3.6 (Table 3, Figure 8) and S1
  (section G.5); figures in the manuscript's terms; nine independent checkers, an editor and two rechecks.
  Every reported sign-flip p is now the exact p over all sign patterns, computed from the stored
  per-location differences (`scripts/g4_confirm_exact_p.py`); the run's Monte Carlo estimates give the same
  Holm decisions in all 25 families. A usage limit stopped the last editor and the response-letter agents;
  the remaining recheck findings were fixed and the point-by-point response written directly, its numbers
  checked against the rendered pages. The revision was merged into the main branch, and every export and
  report figure re-drawn at that clean commit (content identical; only the provenance stamps changed).
- 2026-10-05: round 2 returned: minor revisions (Fable 5.1) and major revisions (GPT-6 Astra), so a third
  round is needed. The round-3 revision (92c19ac and the commits after it): the three school-aged children
  are provisional geometric examples throughout, with the two checks that make them so stated in Section
  2.2 (depth to the MRI boundary; the cortex the source rule leaves out, 3.3-10.0 % in the children against
  7.2 % in the adult, from the new scripts/study_cortex_exclusion.py); the earlier log's "shallow cortex is
  genuine" is narrowed to what the check shows (not a scalp or registration error). The helmet result now
  reads that fitting the helmet to the adult's median gap substantially reduces the smaller heads' extra
  advantage, without isolating the gap as a cause. Also: multi-axis clinical arrays named, the confirmatory
  run's commits and dates and its pilot runs in Section 2.5, full-precision adult detectabilities
  (scripts/export_g2_target_metrics.py), the within-head contrast's aggregation order, a plainer Abstract and
  Highlights, and a point-by-point response to both round-2 reports.
- 2026-10-05: before round 3, an internal pre-check (five independent checkers: response accuracy, coverage of
  every referee point, framing consistency across all pages, numbers, and a mock strict referee; every serious
  finding re-checked by a skeptic) found no wrong number but several overstatements and stale pages. Fixed: the
  fitted helmet now reported with the interaction (+0.40 to +1.10 dB, every interval above zero), which does not
  depend on the adult's reference placement, and with the residual against the adult at top contact; the
  children flagged as provisional wherever their values appear; the templates' scalp convention stated as a
  limitation; the excluded cortex for all nine heads; Figure 7 labelled with its paired estimands and Figure 8
  resized to a journal page; S1, S6, S7, the README, the G2 report and the milestone pages brought in line with
  the manuscript. A second pre-check runs before the package goes to both referees.
- 2026-10-05: second internal pre-check (response accuracy, framing consistency across every page, and two mock
  referees, both predicting minor revisions). Fixed: the reduction stated against each adult reference, the
  templates' scalp convention with its direction, point estimates labelled where intervals span the null, scenario C
  as the most pessimistic assumption rather than a bound, Table 1 with the interaction, the children's skull and
  fiducial rules stated plainly, Figure 8's legend, S1/S6/S7/README wording, and the review labels removed from the
  revision scripts and the metadata of four results (scripts/neutralize_review_labels.py; no number changed).
  Round 3 goes to both referees next.
- 2026-10-05: round 3 (fabd710) returned minor revisions from both referees, the acceptance condition; their wording
  corrections were applied at once (fb7fd2f; 316 tests pass).
- 2026-10-05: final pass on the referees' remaining minor points. Section 3.5 and the Abstract now lead with what fitting the
  helmet restores within each head and what remains against the adult at top contact, each residual tied to its helmet
  construction; the interaction stays in the main text with its per-head table. Section 2.2 is split, with the children's
  skull, fiducials and MRI check as a list and the specific reason their anatomy is in doubt. The MRI check was run on the
  18- and 12-month templates as well (`results/g3b_templates_qc/`, at fb7fd2f): their scalps lie 1.6 and 1.3 mm outside
  their MRI head boundaries, like the 24-month template's, and the 18-month template is classed misregistered (its white
  surface fits its T1 best after a 2.3-mm shift); this is stated in Sections 2.2, 3.5 and 4.3 and S1 D.3, with the two
  scaled adults shown to carry the fitted-helmet result without the templates. Also: the OPM coverage rule and the three
  adult gaps (Section 2.3); the spike detectors' whitener, the oracle's false-positive probability, the S50 rule and the
  exact p (Section 2.5); 7,657 of 7,661 targets; scenario B's factor in the main text; Table 4's Feys et al. (2022)
  analogue on the school-age scaled adult; keywords; the sample dataset's terms quoted with their source; the history
  rewrite moved to Data and code availability; one phrase for the equal-site-count control. Four independent checkers
  (numbers, referee coverage, coherence, code and docs) found no wrong number; their wording findings were fixed. Left to
  the author: funding, competing-interests and contributions statements, an archive DOI, print-resolution figure files,
  and confirming the sample dataset's terms with its providers. Next: deploy, check the live site, and a confirmation
  round with both referees on the deployed version.
- 2026-10-05: f167fd3 deployed (gh-pages af1ef8e; live check: 14 pages, 240 internal assets, 139 anchors, 0 failures;
  39 external links, 31 load and 8 publisher pages refuse scripts: in a browser the MDPI page loads, and the Wiley and
  RSNA DOIs resolve to their article pages, which show those publishers' bot check). Round 4 on that version: minor
  revisions from both referees. Their round-4 corrections applied: the Abstract cut to 249 words (NeuroImage's limit is
  250; Human Brain Mapping sets none) with every condition kept; the 18-month template flagged in Tables 2 and 3 and the
  Figure 6 to 8 captions, with its MRI sections added to Figure S6; the nine anatomies described as not independent
  anatomical replications; the children's skull rule as implemented, with the single-compartment check beside it;
  Section 3.5's first paragraph split; Section 3.6's checks with each array's false-event rates and the equal-rate tests;
  the triaxial control with noisier tangential axes; which way a displaced cortex would move the 18-month template's Δ;
  "would test", not "would decide", for a measured OPM covariance; the exact-p file as the primary source of the
  confirmatory p values on S8; the dataset versions; the ethics sentence under Declarations; and
  docs/literature/model_counterparts.json marked as superseded by Table 4. Left to the author: funding, competing
  interests and contributions, an archive DOI, print-resolution figure files, Table 4's journal layout, and confirming
  the sample dataset's terms with its providers.
- 2026-10-05: two independent checkers verified the round-4 corrections (104b086): every new number matches its source;
  their wording findings were fixed (f60e083): the Abstract's spike result keeps "fixed adult helmet" and the
  localization result "exploratory" (249 words); the 18-month template's flag uses the MRI check's verdict everywhere;
  the independence of the anatomies and their shared assumptions are stated separately; the equal-rate test is called
  approximate and uncorrected, and the rate matching is to the nominal 1 per minute; the superseded comparison record's
  note lists how it differs from Table 4. Next: tests, deploy, live check.
- 2026-10-05: final deployment. acb4752 (318 tests pass) deployed as gh-pages f8292c3 and checked live: 14 pages, 241
  internal assets, 139 anchors, 0 internal failures; 39 external links, 31 load and the same 8 publisher pages refuse
  scripts. Final referee verdicts: round 4, on the deployed f167fd3, minor revisions from Claude Fable 5.1 (maximum
  effort) and from GPT-6 Astra (maximum reasoning, Codex CLI), as in round 3 (fabd710); the live version adds only the
  corrections they requested in round 4, checked by two independent checkers. Open for the author: funding,
  competing-interests and contributions statements, an archive DOI, print-resolution figure files, Table 4's journal
  layout, confirming the MNE sample dataset's terms with its providers, and the reference details at proof.

# NeuroImage manuscript

**Goal (2026-10-05).** Rewrite the published report as a manuscript ready for submission to NeuroImage, without
re-analysis: `paper/` holds the LaTeX source and the compiled PDF (title, highlights, abstract, keywords, Introduction,
Materials and methods, Results, Discussion, Conclusion, references, print-resolution multi-panel figures and tables),
written to a one-page style sheet derived from ten open-access NeuroImage MEG/OPM-MEG papers of 2024-2026 and the
current Guide for Authors; every number matches the analysis outputs. Done when Claude Fable 5.1 (maximum effort) and
GPT-6 Astra (maximum reasoning, Codex CLI), each seeing only the manuscript, the ten papers and the analysis outputs,
return accept or minor revisions in the same round.

| Step | Status | Notes |
|---|---|---|
| 1. Reference papers and Guide for Authors | done | ten NeuroImage papers (PMC Open Data author manuscripts: PDF, XML, text, figures; 102 files, 75.8 MB, checksums verified) and the Guide for Authors (saved by the owner) in `refs/`, which is not committed |
| 2. Style sheet | done | one page (`paper/style_sheet.md`, removed from the release in 8c77a11 and kept in the history), from eleven parallel analyses and a synthesis |
| 3. Figures at print resolution | done | `paper/export_figures.py`: 10 main and 15 supplementary figures as vector PDFs (raster layers at 600 dpi) from the stored results; bold capital panel letters; American spelling; explanatory notes and figure titles moved out of the images into the captions |
| 4. Manuscript and supplementary material | done | `paper/manuscript.tex.j2` and `paper/supplementary.tex.j2` (+ `paper/supp/`), filled from the facts by `paper/build_paper.py`, which also checks highlights (3-5, at most 85 characters) and the abstract (at most 400 words, by the owner's decision; the Guide for Authors asks for 250); 412 facts in the main text, 2,241 in the supplementary material |
| 5. Internal checks | done | five independent checkers (fidelity to the report, style sheet, mock referee, supplementary material, figures and captions); their findings applied; mean sentence 26.0 words |
| 6. Referee rounds | done | round 1 (924c2c6): major revisions from both; round 2 (ea8d5d9): major revisions from both; round 3 (88267e6): major revisions from both; round 4 (16c3e6a): **minor revisions from both** (Claude Fable 5.1 and GPT-6 Astra), the acceptance condition, with requirements A-D achieved by both; final version 323dd1a with their minor corrections |

**Requirements added by the owner (2026-10-05), each to be confirmed by both referees in the same round (D added later the same day):**

| Requirement | Status |
|---|---|
| A. The paper's prose (style, tone, word selection, sentence structure and length, sentence and paragraph composition, flow) matches the four papers in `prose_example/` (not committed) | main text rewritten: prose profiles of the four papers, a measurable prose guide (46 targets), then a section-by-section rewrite; 44 of 46 targets met (the two misses are bounded by keeping every reported number and by not hedging findings); the adult comparison and the spike detection moved into main-text Tables 3 and 5 so that the Results prose is figure- and table-led |
| B. https://khan-laiba.github.io/opm-squid-meg/ shows the HTML version of the paper | done and live: `paper/build_html.py` (the rendered LaTeX with LaTeX's own numbering and BibTeX's reference list, figures rasterized, equations by MathJax, every local link checked), published by `scripts/deploy_pages.py`; first deployed 2026-10-06 from 8cbb043 (gh-pages bea4a64; live check: 2 pages, 31 internal assets, 235 anchors, 0 failures), redeployed with every revision |
| C. The repository is cleaned and restructured as the accompanying repository of a NeuroImage paper | done (8c77a11, merged to main): code, configurations, results, manuscript sources, documentation and tests; the report, its site builder, PLAN.md and other process material removed (kept in the history); CITATION.cff; a map from every figure and table to its scripts and result files in `scripts/README.md`. README rewritten again on 2026-10-06 at the owner's request, for a scientific reader (see the log) |
| D. The abstract reads as well as the abstracts of the four example papers (style, flow, prose); written last, after the whole final paper is read; up to 400 words allowed by the owner (the Guide for Authors asks for 250) | rewritten in full on 2026-10-06 at the owner's request, after rereading the whole paper, in the flow of the abstract of 2026.08.17.744953 (premise; the open question; the factors examined; findings in words with their conditions; "Overall, ..."); 397 words; five new highlights; to be judged by both referees |

## Manuscript log

- 2026-10-05: branch `neuroimage-manuscript`. ScienceDirect refuses automated access, so the ten papers come from the
  PMC Open Data bucket (NIH-deposited author manuscripts of the accepted papers, CC BY/BY-NC/BY-NC-ND): Matsubara 2026,
  Núñez Ponasso 2025, Wartman 2025, Arif 2025, Youssofzadeh 2026, Jiao 2024, Ward 2025, Pulliam 2024, Jia 2025, Wang
  2025. The owner saved the Guide for Authors. Toolchain: tectonic (XeTeX) with Elsevier's elsarticle class and TeX Gyre
  Termes. Title page and declarations as decided by the owner: Laiba Khan (Lexington High School) and Mainak Jas
  (Martinos Center, MGH/HMS; corresponding author); no specific funding; no competing interests; CRediT roles; a
  generative-AI declaration.
- 2026-10-05: manuscript drafted. Main text written to the style sheet from the report (structure 1-5, Discussion with
  an unnumbered summary and the limitations in one paragraph, no bullet lists, no provenance material); the
  supplementary material converted section by section by four parallel writers, each fragment checked to compile
  alone. Sentences split to a mean of 25.8 words. Floats carry no hyperlinks (the PDF driver misplaces them).
- 2026-10-05: internal checks. Five independent checkers read the manuscript, the supplementary material and the
  figures against the reviewed report, the style sheet and the reference papers. Applied: conditions restored to the
  highlights, abstract, Discussion and Conclusions (children provisional, anatomies not independent replications,
  idealized detector, exploratory localization); a causal reading ("helmet fit rather than head size") withdrawn; the
  within-head contrast and the interaction stated as independent of the assumed OPM noise; the break-even OPM noise
  stated as the testable prediction; tables reordered so that Table 2 (estimands) follows the definitions; captions
  shortened; a wider text block so that figures print near their drawn size; the scaled adults and Fig. 7's helmets
  given colors no array uses; Section S8.3 (literature search and parameter provenance) restored from the report.
- 2026-10-05: referee round 1 on 924c2c6. Claude Fable 5.1 (maximum effort): major revisions. GPT-6 Astra (maximum
  reasoning, Codex CLI): major revisions. Both found no mismatched number. Their shared points: lead with what is robust
  (array design, OPM noise, comparator and metric decide the sign; the break-even noise as the testable prediction; the
  helmet-fit interaction, independent of the OPM noise), present the measured-Neuromag scenario with the primary
  result, treat the parcel intervals as descriptive, present the declared spike run as a verification under the model,
  base the pediatric conclusion on the scaled adults, and move the literature table to the supplement. Requests for new
  simulations (a clinical triaxial array, spikes under other noise models or in the fitted helmet, repaired anatomies,
  longer threshold calibration) are outside this rewrite; the claims were narrowed instead and the existing results
  that bear on them brought into the main text. Round 2 revision: new abstract and highlights, ten main figures (the
  stacked figures split; a helmet-fit figure that plots the interaction), Table 2 pruned, Section S9 (Table S19), the
  supplement tightened.
- 2026-10-05: round 2 revision complete and sent to both referees with a point-by-point response. The supplement was
  condensed by four parallel writers (Sections S3, S4 and S7 by about a sixth; repeats of tables and of the main text
  removed; post hoc choices disclosed once, in Section S8.2), with Figs. S3, S5, S14 and S15 enlarged. Every
  supplementary section, figure and table the main text cites was checked for number and content. The build is clean:
  main text 29 pages, supplementary material 56 pages, no undefined references or float errors.
- 2026-10-05: round 2 (ea8d5d9) returned major revisions from both referees. Both again found no wrong number and judged
  most first-round points addressed. Their remaining points: the abstract in the PDF was cut short by an unescaped
  percent sign (the build now rejects a bare % in the body and checks that the abstract's last words reach the PDF);
  the break-even noise levels and the scenario range must be stated as conditional on the modeled SQUID noise; the
  misregistered 18-month template and the three children must be separated from the principal pediatric evidence; the
  interaction must be defined as computed (a paired median of target-level contrasts, not a difference of medians);
  the bracketing of clinical arrays and the no-reversal inference for spikes must go; the detection comparison holds
  at nominal false-event rates; and the Methods need the forward-model, covariance-estimation and estimand details.
  Round 3 revision under way, without new simulations.
- 2026-10-05: round 3 revision, main text. The round-2 points were addressed in wording and presentation only: the
  break-even levels stated as conditional on the modeled Neuromag noise (under the measured noise the dense array's
  level would lie a little above 15 fT/√Hz); the principal smaller heads separated from the four heads of a new Section
  S10; the interaction defined as computed (Eq. 3); new title. Then the prose rewrite of requirement A: Results
  paragraphs led by figures and tables, new Table 3 (the adult comparison under every condition, comparator, measure,
  model variant and scenario) and Table 5 (spike detection by detector and thresholds), the limitations kept in one
  Discussion paragraph, sentences of 22 words on average. Figure corrections (labels, legends, the misregistered
  18-month template marked in Figs. S16-S18) under way. Owner, later the same day: the abstract is to be rewritten last,
  after the whole final paper is read, and both referees gate it (requirement D).
- 2026-10-05/06: two internal pre-checks of the rewritten main text (numbers, referee coverage, a mock referee,
  consistency and layout, prose against the four example papers, the supplement). They found three qualitative claims
  of the rewrite that the stored values contradict (an "only" reversal, "the largest effect of any single change", an
  incomplete list of the site-matched array's advantages), mislabeled Table 3 rows, and conclusions worded more broadly
  than the evidence (the helmet-fit result holds with each head's own, smaller dense array; proximity was never ranked
  against the other factors; under the measured Neuromag noise an equal-SNR depth returns near the sphere's). All were
  corrected in wording, with every number still read from the results; the figures were corrected (labels, legends,
  overlaps, minimum type of 6.5 pt, the misregistered template marked); the supplement is being aligned with the main
  text. The abstract and highlights were then written last. Unit tests: 318 pass.
- 2026-10-06: repository reorganized as the paper's code and data release (8c77a11): the report, its site builder,
  PLAN.md and process notes removed (all kept in the history); new README, CITATION.cff, content licence and a map
  from every figure and table to its scripts and result files (`scripts/README.md`). The HTML version of the paper
  was deployed to GitHub Pages from 8cbb043 (gh-pages bea4a64; live check: 2 pages, 31 internal assets, 235 anchors,
  0 failures). Unit tests: 299 pass (the site builder's tests went with it).
- 2026-10-06: owner: the abstract and the README were not good enough ("Think like a scientist. This is a top
  academic paper."), and "record larger fields" is not scientific language. The whole paper was reread and the
  abstract rewritten from scratch in the flow of the abstract of 2026.08.17.744953: the premise (on the scalp, the
  field of cortical sources is stronger than at the SQUIDs), the open question (clinical comparisons found a higher
  OPM SNR in some studies but not in others), the spherical prediction, the factors examined and the design, the
  findings with their conditions, and an "Overall, ..." conclusion; 396 words, every number a fact. The five
  highlights follow it. The README was rewritten for a scientific reader: an overview of the question and design
  with Figs. 1 and 2, the main results with their conditions and figure and table references, system requirements,
  installation, the external data with their terms, reproduction from the stored results and from the data (with
  measured run times), the declaration of the pre-specified spike run, tests, citation, licence and contact.
  Two independent pre-checks followed (abstract and highlights against the paper, the facts and the example
  abstracts; README against the paper, the repository, the git history and the run logs). They found the
  measured-noise ratio stated without its assumption (Neuromag measured, OPM modeled), "identified" used for a
  detection endpoint, two loose verbs, README run times summed from overlapping timers (the pre-specified spike run
  took about 3 h, not 4 h), result files stamped with pre-rewrite hashes without a pointer to docs/commit_map.tsv,
  and missing data terms; all corrected, together with the same colloquial verbs in the main text and the
  supplement (abstract 397 words). Committed as 88267e6 (main fast-forwarded), deployed as gh-pages 55761a3 and
  checked live (2 pages, 31 internal assets, 235 anchors, 0 failures; the new abstract is served). Unit tests: 299
  pass; facts: 68,944, 0 problems. Round 3 sent to both referees on 88267e6.
- 2026-10-06: round 3 (88267e6) returned major revisions from both referees; neither found a wrong number (Fable 5.1
  checked 26 groups of values, GPT-6 Astra 26). Requirements: abstract achieved (both); HTML achieved (Fable 5.1),
  not achieved (GPT-6 Astra: Table S19's continuation numbered S20 and its notes lost); repository achieved (Fable
  5.1), not achieved (GPT-6 Astra: three superseded v1-array diagnostics still in results/); prose achieved (GPT-6
  Astra), not achieved (Fable 5.1: uniformly medium sentences, enumerative Results, coined labels). Both: the
  scenario-D break-even level was inferred, not computed. Fable 5.1 also asked what "head size" means for scaled heads.
- 2026-10-06: round 3 revision (b85dcdf, 16c3e6a), no new simulations. The scenario-D break-even values were removed
  (the stored results support only the direction at 15 fT/√Hz); the head-size caveat added to Section 4.2 and the
  Limitations; the HTML converter fixed for continued tables, grouped notes and line breaks in cells (every table of
  both documents checked against the PDF; three new unit tests); the v1 diagnostics and their loading removed; Figs. 4D
  and 5A relabeled; about 40 targeted corrections. The Methods and Results were then revised for prose by six
  writer-verifier pairs, block by block, with every fact, float and reference locked; two independent pre-checks
  followed (a claims and consistency audit, which found one misattributed quantity in the new head-size paragraph, and
  a prose judge applying Fable 5.1's criteria, which judged the requirement achieved after trimming re-definitions in
  the Results). All findings applied. Deployed as gh-pages 0e13a01 and checked live (0 failures). Unit tests: 302
  pass; facts: 68,864, 0 problems. Round 4 sent to both referees on 16c3e6a.
- 2026-10-06: round 4 (16c3e6a) returned **minor revisions from both referees** (Claude Fable 5.1 and GPT-6 Astra), the
  acceptance condition of the goal. Both judged all four requirements achieved (prose against the example papers, the
  HTML version, the accompanying repository, the abstract), found no remaining major issue and again no wrong number.
  Their minor points were applied in 323dd1a: the measured-covariance interpretation stated exactly, the known noise
  covariance called the true covariance throughout, the scenario-D break-even direction qualified as a point estimate
  in the Discussion, the abstract's spike range tied to the adult and the smaller heads in the adult helmet, notes on
  the rounded bounds of Table 4 and the pooled row of Table 5, and small repository cleanups; an independent verifier
  checked every correction against both reports and the result files. Not applied, by the goal's or the authors'
  choice: splitting the Limitations paragraph (the goal keeps it as one), the full text of Ren et al. (2025) (not
  readable here; needs institutional access), consolidating supplementary tables. Unit tests: 302 pass; facts:
  68,864, 0 problems. Deployed as gh-pages 9c798ca.
