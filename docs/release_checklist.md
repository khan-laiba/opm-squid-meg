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

Last checked 2026-10-02: `private: true`, `visibility: private`, `has_pages: false`; the Pages
endpoint returns 404 (no site). Repository privacy does not establish website privacy: check Pages
separately after any change. GitHub's documentation (read 2026-10-02, DOC-PAGES in the register)
allows access-controlled Pages sites only on GitHub Enterprise Cloud; this repository belongs to a
user account, so a Pages site, if ever enabled, would be public.

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
- Not committed: reference PDFs, MNE sample data and anatomy, the infant templates, the school-aged
  children (OpenNeuro ds005234), fsaverage,
  caches (`.gitignore`; `git ls-files` lists only generated figures as PDFs).
- No secrets: `git grep -I -i -E "api[_-]?key|secret|token|passw"` finds only the goal text and
  this checklist.
- The local report builds and link-checks (`scripts/build_site.py`, every link and download
  checked; `tests/test_site.py`).

## 3. Owner decisions before any public release

1. **Licence.** The repository has no licence file; without one the code is not reusable by
   others. Choose a licence (and whether results and figures carry a different one).
2. **Identifying metadata and history.** All commits carry the author name and institutional e-mail
   address, and the history (not the current tree) also holds two older text blobs with an e-mail
   address and the earlier PDFs that embed the licensed artwork fonts (e.g. `git show
   81168f3:legacy/figure3_output/Figure3_replicated.pdf | pdffonts -`). Publishing the history as it
   is publishes all of these, so a public release needs a fresh export or a rewritten history
   (which needs a force-push); the choice is the owner's. No force-push or history rewrite has been
   done. This is a release blocker, not a defect of the private results.
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
4. **School-aged children's derivatives** (OpenNeuro ds005234, Fadeev et al. 2024): `results/g3b/g3b_targets_child*.csv`,
   `results/g3b/school_subjects_preparation.json`, `results/g3b/school_anatomy_checks.json`, the G3B summary, report
   and figures that include children A-C, `results/g4/*child*`, the pediatric G4 cross-reading summaries and reports
   (`g4_pediatric_comparison.json`, `G4_pediatric_report.md`, `g4_matched_rate.json`, `G4_matched_rate_report.md`,
   `g4_fit_failures.json`, `G4_fit_failures_report.md`) and the report pages built from them. The dataset's metadata
   say CC0, its acknowledgement text CC BY; cite it and confirm. `configs/school_subjects_manifest.json` lists the
   source files only (public S3 object versions, sizes, SHA-256; no data).
5. **fsaverage derivatives** (`legacy/realistic_head_output/`): FreeSurfer's fsaverage; confirm
   the terms for derived figures and per-vertex tables.
6. **Published figure raster.** `legacy/figure3_output/Figure3_published_extracted.png`,
   `Figure3_overlay.png`, `Figure3_side_by_side.png` and `Figure3_difference.png` contain the
   published Jas et al. Fig. 3 raster (bioRxiv, doi 10.64898/2026.08.17.744953; the preprint's own
   header states a CC-BY 4.0 International licence, checked 2026-10-01). Keep the attribution with
   them, or leave them out of the release.
7. **Digitised figure data.** `scripts/digitise_hunold_fig6.py` and its output reproduce a
   waveform digitised from Hunold et al. (2016) Fig. 6; confirm that publishing the digitised
   values is acceptable.
8. **Local references.** `docs/audit.md` and `goal_condition.txt` name the local working folders,
   and `legacy/realistic_head_output/realistic_head_results.npz` stores the local SUBJECTS_DIR of
   the original run; edit or drop before publication if the folder names should not appear.
9. **Literature extractions.** `docs/literature/{jas2026,hunold2016,goldenholz2009}.md` are long
   paraphrased extractions, with parameter tables and digitised figure values, of three papers, two
   of them paywalled; confirm that publishing them is acceptable, or keep them private and cite page
   numbers only.
10. **The goal text.** `GOAL.md` and its copy `OPM_SQUID_Adult_to_Pediatric_Goal.md` name a third
   person and describe a mentoring relationship; decide whether they belong in a public release.
11. **Historic diagnoses.** `results/g2/bem_skin_refinement_v1_arrays.json` and
    `results/g2/near_mesh_check_v1_arrays.json` record a dirty commit (historic v1 diagnoses, labelled
    in `docs/methods.md` and in the report's download list); keep them labelled or drop them.
12. **Visibility and Pages.** Only after 1-11: change visibility and, separately, enable Pages;
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
