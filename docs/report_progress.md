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
| 4. Referee rounds (both referees, same round) | done | round 1: major revisions from both; round 2: minor revisions (Fable 5.1), major revisions (GPT-6 Astra); round 3: minor revisions from both; their minor points applied, and the final version sent to both for confirmation (round 4) |
| 5. Deploy to GitHub Pages and check the live site | in progress | the final version is deployed from the commit that applies the round-3 minor points; live check follows |

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
  points were then applied in a final pass (log, 2026-10-05), and the final version goes to both referees once more
  as a confirmation round (round 4).

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
