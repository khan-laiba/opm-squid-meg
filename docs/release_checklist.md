# Release checklist (prepared, not executed)

Status on 2026-10-01: **release-ready: no** (the owner decisions of section 3 are open);
**publicly deployed: no** (private repository, no Pages site).

The repository stays private and no website is deployed. Making the repository public, enabling
GitHub Pages or publishing the report needs the owner's explicit approval (GOAL.md, G5); release-ready
and publicly deployed are separate statuses. This list records what has been prepared and what
remains the owner's decision. It does not claim that any author of the reproduced papers has
reviewed or approved the study.

## 1. Visibility (verify before and after any change)

```bash
gh api repos/khan-laiba/opm-squid-meg --jq '{private: .private, visibility: .visibility, has_pages: .has_pages}'
```

```bash
gh api repos/khan-laiba/opm-squid-meg/pages
```

Last checked 2026-10-01: `private: true`, `visibility: private`, `has_pages: false`; the Pages
endpoint returns 404 (no site). Repository privacy does not establish website privacy: check Pages
separately after any change.

## 2. Done in preparation

- Figures use open fonts (DejaVu Sans, STIX; shipped with matplotlib). The Fig. 3-style figures
  (`results/g1a/fig3/`, `legacy/figure3_output/`, `legacy/extended_sources_output/`,
  `legacy/realistic_head_output/`) were regenerated with them; `pdffonts` on every tracked PDF
  lists only DejaVu and STIX faces. The licensed artwork fonts (Myriad Pro, Times New Roman) are
  used only with `--artwork-fonts`, for the pixel verification of the Fig. 3 replica, and such
  outputs are not committed (D-REL-FONTS). With them the replica passes 7/7 checks and the
  negative controls (`legacy/figure3_output/verification_report_artwork_fonts.json`); with open
  fonts the four geometric checks pass and the three typography checks do not
  (`verification_report.json`).
- Output summaries record file names relative to the summary instead of local absolute paths
  (legacy scripts), and the realistic-head summary no longer records the local SUBJECTS_DIR.
- Not committed: reference PDFs, MNE sample data and anatomy, the infant templates, fsaverage,
  caches (`.gitignore`; `git ls-files` lists only generated figures as PDFs).
- No secrets: `git grep -I -i -E "api[_-]?key|secret|token|passw"` finds only the goal text and
  this checklist.
- The local report builds and link-checks (`scripts/build_site.py`, every link and download
  checked; `tests/test_site.py`).

## 3. Owner decisions before any public release

1. **Licence.** The repository has no licence file; without one the code is not reusable by
   others. Choose a licence (and whether results and figures carry a different one).
2. **Identifying metadata.** All commits carry the author name and institutional e-mail address
   (no tracked file spells it out). Publishing the history as it is publishes them;
   the alternatives (a fresh export, or a rewritten history, which needs a force-push) are the
   owner's call. No force-push or history rewrite has been done.
3. **Infant-template derivatives.** Figures and tables derived from the O'Reilly et al. (2021)
   templates (LGPL-2.1 repository; built from the Neurodevelopmental MRI Database, Richards et al.
   2016, which has its own terms): every `results/g3b` figure and target table that shows a
   template (`Figure_G3B_maps_infant*.png`, `Figure_G3B_usefulness_infant*.png`,
   `Figure_G3B_geometry.png`, `Figure_G3B_depth.png`, `Figure_G3B_delta.png`,
   `Figure_G3B_placements.png`, `g3b_targets_infant*.csv`), the G3B summary and report
   (`g3b_summary.json`, `G3B_report.md`), `results/g4/*infant*`, `results/g4/G4_pediatric_report.md`,
   `results/g4/g4_pediatric_comparison.json`, the motion results (`g4_motion_summary.json`,
   `G4_motion_report.md`, `g4_motion_timecourse_example.csv`, `Figure_G4_motion.png`), and the
   report pages built from them (pediatric, epilepsy, G3B and motion reports, downloads). Confirm
   redistribution; both papers are cited wherever these appear.
4. **fsaverage derivatives** (`legacy/realistic_head_output/`): FreeSurfer's fsaverage; confirm
   the terms for derived figures and per-vertex tables.
5. **Published figure raster.** `legacy/figure3_output/Figure3_published_extracted.png`,
   `Figure3_overlay.png`, `Figure3_side_by_side.png` and `Figure3_difference.png` contain the
   published Jas et al. Fig. 3 raster (bioRxiv 2026.08.17.744953; CC-BY 4.0 as recorded in the
   README). Keep the attribution with them, or leave them out of the release.
6. **Digitised figure data.** `scripts/digitise_hunold_fig6.py` and its output reproduce a
   waveform digitised from Hunold et al. (2016) Fig. 6; confirm that publishing the digitised
   values is acceptable.
7. **Local references.** `docs/audit.md` and `goal_condition.txt` name the local working folders,
   and `legacy/realistic_head_output/realistic_head_results.npz` stores the local SUBJECTS_DIR of
   the original run; edit or drop before publication if the folder names should not appear.
8. **Visibility and Pages.** Only after 1-7: change visibility and, separately, enable Pages;
   re-run section 1.

## 4. Re-run before release

```bash
.venv/bin/python -m unittest discover -s tests -t .
```

```bash
.venv/bin/python scripts/build_site.py
```

```bash
for f in $(git ls-files '*.pdf'); do pdffonts "$f" | tail -n +3 | awk '{print $1}'; done | sed 's/^[A-Z]*+//' | sort -u
```
