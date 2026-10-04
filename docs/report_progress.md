# Report progress

**Goal.** Publish the pediatric epilepsy research report as a GitHub Pages site built from the
analysis already in this repository, as one coherent story that starts from Jas et al. (2026) and
builds step by step to pediatric epilepsy, following the story plan and the supervisor's guidance
(the brain regions covered by OPM epilepsy studies, OPM vs SQUID by region across head sizes, the
number of sensors). Every number and figure must trace to an analysis output in the repository.

**Done when** the site builds and deploys, the live URL serves the final report with every link
and figure resolving, and two independent referees, Claude Fable 5.1 (maximum effort, a fresh
subagent) and GPT-6 Astra (maximum reasoning, Codex CLI), each return "accept" or "minor
revisions" in the same round. Each referee sees only the files.

**Site:** https://khan-laiba.github.io/opm-squid-meg/ (Pages enabled 2026-10-04; nothing deployed
yet). The repository stays private; the Pages site is public.

## Status

| Step | Status | Notes |
|---|---|---|
| 1. Foundations: number provenance, new figures from stored outputs, verified literature table, site integration | done | 31,150 sourced facts; figures R1-R10; 15 studies (14 verified in full text); 197 tests pass |
| 2. Draft the report as one narrative | in progress | one writer, then nine checkers and one editor |
| 3. Internal verification: numbers, claims against evidence, methods completeness, build, tests, links | pending | |
| 4. Referee rounds (both referees, same round) | pending | |
| 5. Deploy to GitHub Pages and check the live site | pending | |

## Referee rounds

None yet.

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
