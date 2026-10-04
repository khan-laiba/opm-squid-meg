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
| 1. Foundations: number provenance, new figures from stored outputs, verified literature table, site integration | in progress | parallel agents |
| 2. Draft the report as one narrative | pending | |
| 3. Internal verification: numbers, claims against evidence, methods completeness, build, tests, links | pending | |
| 4. Referee rounds (both referees, same round) | pending | |
| 5. Deploy to GitHub Pages and check the live site | pending | |

## Referee rounds

None yet.

## Log

- 2026-10-04: Pages enabled on the repository (served publicly, repository private); Codex CLI
  installed and logged in; the publishable files checked for local paths and e-mail addresses (none);
  shared interfaces added (`scripts/report_facts.py`, `scripts/report_style.py`).
