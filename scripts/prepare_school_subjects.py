#!/usr/bin/env python3
"""Prepare the school-aged children of OpenNeuro ds005234 for G3B and the pediatric G4 studies.

Input (git-ignored, data/external/school_subjects/, see its MANIFEST.json): for each child its own
FreeSurfer white and sphere surfaces, aparc annotations and dense MRI scalp (lh.seghead), and its
watershed BEM surfaces, which the dataset stores in another subject's folder (configs/
g3b_pediatric.toml, [[anatomy.school]] bem_folder; checked here by the files' volume information:
the same source volume and c_ras as the child's own surfaces).
Output, laid out as MNE packages the infant templates (data/external/school_subjects/<subject>/bem/):
  <subject>-5120-5120-5120-bem.fif  3-layer BEM, 5,120 triangles per surface: the watershed outer
                                    skin as the head surface (put on the scalp at loading,
                                    A-BEM-CONFORM) and a modelled skull (A-BEM-CHILD: the watershed
                                    inner skull lies just below the scalp)
  <subject>-head.fif                the dense MRI scalp
  <subject>-oct-6-src.fif           oct-6 surface source space on the white surface (with patch
                                    statistics; the full white surface is the full-resolution cortex)
  <subject>-fiducials.fif           LPA, nasion, RPA (MRI frame) transferred from the adult (A-G3-FID)
and results/g3b/school_subjects_preparation.json (checks; derived quantities only).
Nibabel is not needed (opmsquid.fsio).
"""
from __future__ import annotations

import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import mne  # noqa: E402
import numpy as np  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402

from opmsquid import anatomy, fsio, io, paths  # noqa: E402

OUT = ROOT / "results" / "g3b"
BEM_FILES = ("inner_skull.surf", "outer_skull.surf", "outer_skin.surf")


def adult_reference():
    """The adult's cortex (white vertices inside its inner skull) and digitised fiducials (MRI frame)."""
    a = anatomy.load_sample()
    cortex = anatomy.full_resolution(a)
    info = mne.io.read_info(paths.SAMPLE_MEG / "sample_audvis_raw.fif", verbose=False)
    names = {FIFF.FIFFV_POINT_LPA: "lpa", FIFF.FIFFV_POINT_NASION: "nasion", FIFF.FIFFV_POINT_RPA: "rpa"}
    fids = {names[d["ident"]]: mne.transforms.apply_trans(a.trans, d["r"]) for d in info["dig"]
            if d["kind"] == FIFF.FIFFV_POINT_CARDINAL}
    return cortex.rr[cortex.valid], fids


def prepare(child: dict, cfg: dict, adult_cortex: np.ndarray, adult_fids: dict) -> dict:
    root = paths.EXTERNAL / anatomy.SCHOOL_SUBJECTS
    name = child["subject"]
    print(f"{child['key']} ({name}) ...", flush=True)
    own = root / name / "surf"
    bem_src = root / child["bem_folder"] / "surf"
    white = [fsio.read_geometry(own / f"{h}.white", read_metadata=True) for h in ("lh", "rh")]
    sphere = [fsio.read_geometry(own / f"{h}.sphere") for h in ("lh", "rh")]
    seg_rr, seg_tris, seg_info = fsio.read_geometry(own / "lh.seghead", read_metadata=True)
    bem = {f: fsio.read_geometry(bem_src / f, read_metadata=True) for f in BEM_FILES}
    cras = np.asarray(white[0][2]["cras"])
    for f, (_, _, info) in [("rh.white", white[1]), ("lh.seghead", (None, None, seg_info))] + list(bem.items()):
        if not np.allclose(info["cras"], cras, atol=1e-3):
            raise ValueError(f"{name}: {f} has c_ras {info['cras']}, the white surface {cras}")
    for f, (_, _, info) in bem.items():
        if f"/{name}/" not in info["filename"]:
            raise ValueError(f"{name}: {f} from {child['bem_folder']} was made from {info['filename']}")
    scalp = anatomy._outward(mne.surface.complete_surface_info(dict(rr=seg_rr / 1e3, tris=seg_tris, np=len(seg_rr),
                                                                    ntri=len(seg_tris)), copy=False, verbose=False))
    cortex_rr = np.vstack([w[0] for w in white]) / 1e3
    ico4 = {f: mne.bem._ico_downsample(dict(rr=rr, tris=tris), 4) for f, (rr, tris, _) in bem.items()}  # [mm], 5,120 triangles
    inner_ws = mne.surface.complete_surface_info(dict(id=FIFF.FIFFV_BEM_SURF_ID_BRAIN, coord_frame=FIFF.FIFFV_COORD_MRI,
                                                      rr=ico4["inner_skull.surf"]["rr"] / 1e3, tris=ico4["inner_skull.surf"]["tris"],
                                                      np=len(ico4["inner_skull.surf"]["rr"]), ntri=len(ico4["inner_skull.surf"]["tris"])),
                                                 copy=False, verbose=False)
    to_scalp = anatomy.MeshDistance(dict(rr=scalp.rr, tris=scalp.tris))
    ws_depth = to_scalp.unsigned(inner_ws["rr"]) * 1e3
    inner, outer, skull = anatomy.model_skull(inner_ws, scalp, cortex_rr, depth=cfg["anatomy"]["child_skull_depth_mm"] * 1e-3,
                                              min_cortex=cfg["anatomy"]["child_skull_min_cortex_mm"] * 1e-3)
    skin_ws = mne.surface.complete_surface_info(dict(id=FIFF.FIFFV_BEM_SURF_ID_HEAD, coord_frame=FIFF.FIFFV_COORD_MRI,
                                                     rr=np.array(ico4["outer_skin.surf"]["rr"], float) / 1e3,
                                                     tris=ico4["outer_skin.surf"]["tris"], np=len(ico4["outer_skin.surf"]["rr"]),
                                                     ntri=len(ico4["outer_skin.surf"]["tris"])), copy=False, verbose=False)
    (skin, _, _), conform = anatomy.head_on_scalp([skin_ws, outer, inner], scalp)  # A-BEM-CONFORM before writing
    skin = dict(rr=np.asarray(skin["rr"], float) * 1e3, tris=skin["tris"])
    surfs = mne.bem._surfaces_to_bem(
        [dict(rr=inner["rr"] * 1e3, tris=inner["tris"]), dict(rr=outer["rr"] * 1e3, tris=outer["tris"]),
         dict(rr=np.array(skin["rr"], float), tris=skin["tris"])],
        [FIFF.FIFFV_BEM_SURF_ID_BRAIN, FIFF.FIFFV_BEM_SURF_ID_SKULL, FIFF.FIFFV_BEM_SURF_ID_HEAD], [0.3, 0.006, 0.3], ico=None)
    head = mne.bem._surfaces_to_bem([dict(rr=seg_rr, tris=seg_tris)], [FIFF.FIFFV_BEM_SURF_ID_HEAD], [1.0], incomplete="warn")[0]
    src = fsio.surface_source_space([w[:2] for w in white], sphere, name, "oct6")
    fids, fid_report = anatomy.transfer_fiducials(adult_cortex, adult_fids, cortex_rr, scalp.rr)
    bem_dir = root / name / "bem"
    bem_dir.mkdir(parents=True, exist_ok=True)
    mne.write_bem_surfaces(bem_dir / f"{name}-5120-5120-5120-bem.fif", surfs, overwrite=True, verbose=False)
    mne.write_bem_surfaces(bem_dir / f"{name}-head.fif", head, overwrite=True, verbose=False)
    mne.write_source_spaces(bem_dir / f"{name}-oct-6-src.fif", src, overwrite=True, verbose=False)
    ident = {"lpa": FIFF.FIFFV_POINT_LPA, "nasion": FIFF.FIFFV_POINT_NASION, "rpa": FIFF.FIFFV_POINT_RPA}
    mne.io.write_fiducials(bem_dir / f"{name}-fiducials.fif",
                           [dict(kind=FIFF.FIFFV_POINT_CARDINAL, ident=ident[k], r=np.asarray(v, float)) for k, v in fids.items()],
                           coord_frame=FIFF.FIFFV_COORD_MRI, overwrite=True, verbose=False)
    sub = anatomy.load_school(name)  # the loader's own checks (A-BEM-CONFORM)
    cor = anatomy.full_resolution(sub)
    return dict(key=child["key"], subject=name, age_y=child["age_y"], bem_folder=child["bem_folder"],
                n_white_vertices=[int(len(w[0])) for w in white], n_scalp_vertices=int(len(seg_rr)),
                watershed_inner_skull_scalp_depth_mm_p10_p50_p90=[float(x) for x in np.percentile(ws_depth, [10, 50, 90])],
                skull_model=skull, fiducials=fid_report, head_conform=conform, head_conform_at_loading=sub.head_conform,
                cortex_inside_inner_skull_share=float(cor.valid.mean()), cortex_usable_share=float(cor.usable.mean()),
                oct6_sources=[int(s["nuse"]) for s in src])


def main():
    cfg = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    mne.set_log_level("WARNING")
    adult_cortex, adult_fids = adult_reference()
    out = dict(status="NEW (G3B: school-aged children of OpenNeuro ds005234 prepared; A-BEM-CHILD, A-G3-FID, A-BEM-CONFORM)",
               children=[])
    for child in cfg["anatomy"]["school"]:
        rep = prepare(child, cfg, adult_cortex, adult_fids)
        out["children"].append(rep)
        print(f"{child['key']} ({child['subject']}): skull moved {rep['skull_model']['moved_median_mm']:.1f} mm (median), "
              f"inner skull depth {rep['skull_model']['scalp_depth_after_mm'][1]:.1f} mm, cortex usable {rep['cortex_usable_share']:.3f}, "
              f"fiducial fit {rep['fiducials']['cortex_fit_median_mm']:.1f} mm", flush=True)
    io.write_json(out, OUT / "school_subjects_preparation.json")


if __name__ == "__main__":
    main()
