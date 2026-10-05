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

**Site:** https://khan-laiba.github.io/opm-squid-meg/ (served from the `gh-pages` branch; the version
under review is live, and the accepted version will replace it). The repository stays private; the
Pages site is public.

## Status

| Step | Status | Notes |
|---|---|---|
| 1. Foundations: number provenance, new figures from stored outputs, verified literature table, site integration | done | 31,150 sourced facts; figures R1-R10; 15 studies (14 verified in full text); 197 tests pass |
| 2. Draft the report as one narrative | done | about 8,500 words of main text, 10 figures, 3 tables, 39 references |
| 3. Internal verification: numbers, claims against evidence, methods completeness, build, tests, links | done | 166 checker issues resolved; supplementary pages cleaned for publication; 197 tests pass |
| 4. Referee rounds (both referees, same round) | in progress | round 1 |
| 5. Deploy to GitHub Pages and check the live site | first deployment done | the round-1 version is live; redeploy after acceptance |

## Referee rounds

- Round 1 (commit f59f76e): both referees reviewing.

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
