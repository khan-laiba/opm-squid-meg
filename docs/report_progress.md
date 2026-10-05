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

**Site:** https://khan-laiba.github.io/opm-squid-meg/ (served from the `gh-pages` branch; the round-1
version is live, and the revised version will replace it). The repository is public (code MIT; report,
figures and result files CC BY 4.0).

## Status

| Step | Status | Notes |
|---|---|---|
| 1. Foundations: number provenance, new figures from stored outputs, verified literature table, site integration | done | 31,150 sourced facts; figures R1-R10; 15 studies (14 verified in full text); 197 tests pass |
| 2. Draft the report as one narrative | done | about 8,500 words of main text, 10 figures, 3 tables, 39 references |
| 3. Internal verification: numbers, claims against evidence, methods completeness, build, tests, links | done | 166 checker issues resolved; supplementary pages cleaned for publication; 197 tests pass |
| 4. Referee rounds (both referees, same round) | in progress | round 1: major revisions from both; revision done; round 2 under review |
| 5. Deploy to GitHub Pages and check the live site | first deployment done | the round-1 version is live; redeploy after acceptance |

## Referee rounds

- Round 1 (commit f59f76e): **major revisions** from both referees. They verified the numbers (no
  number contradicted its result file) and asked for: evidence on the noise model that decides the
  comparison (validation against measured Neuromag noise; near-skull cortex, coloured OPM noise and
  cardiac/ocular fields as sensitivity analyses); a quality check of the school-aged children's MRI
  surfaces; a helmet fitted to each head at the adult's gap, with placement uncertainty; a confirmatory
  spike run with the endpoint fixed in advance; public code; and a shorter manuscript in its own terms.
- Round 2 (commit after the revision's clean re-draw): both referees reviewing the revised manuscript, with the round-1 reports and the point-by-point response.

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
