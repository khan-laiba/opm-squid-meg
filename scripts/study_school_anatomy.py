#!/usr/bin/env python3
"""Checks behind the school-aged children's anatomy (G3B; D-G3-ANAT, A-BEM-CHILD, A-G3-FID).

1. Inner skull below the scalp over the upper head (inner-skull vertices above their mean height in
   the MRI frame, as in scripts/study_child_bem.py; distance to the subject's MRI scalp): the
   adult's and the infant templates' segmented inner skulls, the children's watershed inner skulls
   (5,120 triangles, as used) and their modelled ones (A-BEM-CHILD).
2. A-G3-FID on the templates, whose own fiducials are known: the adult's digitised fiducials
   transferred by the cortex fit (``anatomy.transfer_fiducials``, as for the children) against
   the template's own (distance per fiducial; rotation between the two head frames); and the same
   transfer with the fit made between the scalps instead of the cortices.
3. The MNI alternative on the adult: fsaverage's fiducials mapped through the adult's own
   talairach.xfm (``fsio.mni_fiducials``, as ``mne.coreg.get_mni_fiducials``) against its digitised
   ones.
4. The children's talairach.xfm files: the snapshot's file tree shifts them like the BEMs (the S3
   object behind each file, configs/school_subjects_manifest.json, names its owner). White-surface
   centroid in MNI305 against the adult's (through each one's talairach.xfm), and the median
   distance of the child's MNI cortex to the adult's, for each child with the file stored in its
   own folder and with its own file where it was downloaded (child A's, in sub-Z209's folder).
5. The children's thinnest layers: white surface to MRI scalp, modelled inner to outer skull, and
   outer skull to the BEM head surface (minimum distances).

Output: results/g3b/school_anatomy_checks.json (derived quantities only).
"""
from __future__ import annotations

import json
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import mne  # noqa: E402
import numpy as np  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

import prepare_school_subjects as prep  # noqa: E402
from opmsquid import anatomy, fsio, io, paths  # noqa: E402

OUT = ROOT / "results" / "g3b"
FID_KEYS = {FIFF.FIFFV_POINT_LPA: "lpa", FIFF.FIFFV_POINT_NASION: "nasion", FIFF.FIFFV_POINT_RPA: "rpa"}


def upper_depth(inner_rr: np.ndarray, scalp_rr: np.ndarray, scalp_tris: np.ndarray) -> dict:
    """Median (10th, 90th percentile) distance [mm] of the inner skull's upper vertices to the scalp."""
    rr = np.asarray(inner_rr, float)
    up = rr[:, 2] > rr[:, 2].mean()
    d = anatomy.MeshDistance(dict(rr=np.asarray(scalp_rr, float), tris=np.asarray(scalp_tris))).unsigned(rr[up]) * 1e3
    return dict(n_vertices=int(up.sum()), p10_p50_p90_mm=[float(x) for x in np.percentile(d, [10, 50, 90])])


def brain(subject) -> dict:
    return next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)


def frame_difference(a: dict, b: dict) -> dict:
    """Distance per fiducial [mm] and the rotation [deg] between the head frames that a and b define."""
    def frame(f):
        return mne.transforms.get_ras_to_neuromag_trans(f["nasion"], f["lpa"], f["rpa"])

    r = frame(a)[:3, :3] @ frame(b)[:3, :3].T
    return dict(distance_mm={k: float(np.linalg.norm(np.asarray(a[k]) - np.asarray(b[k])) * 1e3) for k in ("lpa", "nasion", "rpa")},
                head_frame_rotation_deg=float(np.degrees(np.arccos(np.clip((np.trace(r) - 1) / 2, -1, 1)))))


def mni_cortex(white_files) -> np.ndarray:
    """White-surface vertices [mm] in scanner RAS of their own volume (MNI305 after a talairach.xfm)."""
    pts = []
    for f in white_files:
        rr, _, info = fsio.read_geometry(f, read_metadata=True)
        m = fsio.surface_ras_to_scanner(info)
        pts.append(rr @ m[:3, :3].T + m[:3, 3])
    return np.vstack(pts)


def to_mni(scanner_rr: np.ndarray, tal: np.ndarray) -> np.ndarray:
    return scanner_rr @ tal[:3, :3].T + tal[:3, 3]


def main():
    cfg = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    manifest = json.loads((ROOT / "configs" / "school_subjects_manifest.json").read_text())
    mne.set_log_level("WARNING")
    root = paths.EXTERNAL / anatomy.SCHOOL_SUBJECTS
    out = dict(status="NEW check: the school-aged children's skull, fiducials and talairach files (D-G3-ANAT, A-BEM-CHILD, A-G3-FID)",
               upper_head_definition="inner-skull vertices above their mean height (MRI frame z), as in scripts/study_child_bem.py")

    # 1. inner skull below the scalp over the upper head
    depth, thin = {}, {}
    adult = anatomy.load_sample()
    depth["adult"] = upper_depth(brain(adult)["rr"], adult.scalp.rr, adult.scalp.tris)
    for t in cfg["anatomy"]["templates"]:
        s = anatomy.load_template(t)
        depth[t] = upper_depth(brain(s)["rr"], s.scalp.rr, s.scalp.tris)
    for child in cfg["anatomy"]["school"]:
        name = child["subject"]
        rr, tris = fsio.read_geometry(root / child["bem_folder"] / "surf" / "inner_skull.surf")
        ico4 = mne.bem._ico_downsample(dict(rr=rr, tris=tris), 4)
        s = anatomy.load_school(name)
        depth[f"{child['key']}/watershed"] = upper_depth(np.asarray(ico4["rr"], float) / 1e3, s.scalp.rr, s.scalp.tris)
        depth[f"{child['key']}/modelled"] = upper_depth(brain(s)["rr"], s.scalp.rr, s.scalp.tris)
        whites = np.vstack([fsio.read_geometry(root / name / "surf" / f"{h}.white")[0] for h in ("lh", "rh")]) / 1e3
        surf = {b["id"]: b for b in s.bem_surfaces}
        brain_rr, skull_rr = surf[FIFF.FIFFV_BEM_SURF_ID_BRAIN]["rr"], surf[FIFF.FIFFV_BEM_SURF_ID_SKULL]["rr"]
        thin[child["key"]] = dict(
            white_to_scalp_min_mm=float(anatomy.MeshDistance(dict(rr=s.scalp.rr, tris=s.scalp.tris)).unsigned(whites).min() * 1e3),
            inner_to_outer_skull_min_mm=float(anatomy.MeshDistance(surf[FIFF.FIFFV_BEM_SURF_ID_SKULL]).unsigned(brain_rr).min() * 1e3),
            outer_skull_to_head_min_mm=float(anatomy.MeshDistance(surf[FIFF.FIFFV_BEM_SURF_ID_HEAD]).unsigned(skull_rr).min() * 1e3))
    out["inner_skull_scalp_depth_upper_head"] = depth
    out["children_thinnest_layers"] = thin

    # 2. A-G3-FID on the templates
    adult_cortex, adult_fids = prep.adult_reference()
    fid = {}
    for t in cfg["anatomy"]["templates"]:
        bem_dir = paths.EXTERNAL / anatomy.INFANT_SUBJECTS / t / "bem"
        fids, _ = mne.io.read_fiducials(bem_dir / f"{t}-fiducials.fif", verbose=False)
        own = {FID_KEYS[f["ident"]]: np.asarray(f["r"], float) for f in fids if f["ident"] in FID_KEYS}
        s = anatomy.load_template(t)
        cortex = np.vstack([h["rr"] for h in s.src])  # every white-surface vertex, as for the children
        moved, rep = anatomy.transfer_fiducials(adult_cortex, adult_fids, cortex, s.scalp.rr)
        by_scalp, _ = anatomy.transfer_fiducials(adult.scalp.rr, adult_fids, s.scalp.rr, s.scalp.rr)
        fid[t] = dict(frame_difference(moved, own), cortex_fit_median_mm=rep["cortex_fit_median_mm"], scale=rep["scale"],
                      scalp_fit=frame_difference(by_scalp, own))
    out["fiducial_transfer_on_templates"] = fid

    # 3. fiducials from MNI coordinates on the adult
    sd = paths.SUBJECTS_DIR / "sample"
    _, _, info = fsio.read_geometry(sd / "surf" / "lh.white", read_metadata=True)
    adult_tal = fsio.read_talairach(sd / "mri" / "transforms" / "talairach.xfm")
    out["mni_fiducials_on_adult"] = frame_difference(fsio.mni_fiducials(info, adult_tal), adult_fids)

    # 4. the children's talairach.xfm files
    owner = {e["dest"]: re.search(r"/freesurfer/(sub-[^/]+)/", e["source"]).group(1) for e in manifest["files"]}
    tal_files = {d: owner[d] for d in owner if d.endswith("mri/transforms/talairach.xfm")}
    adult_mni = to_mni(mni_cortex([sd / "surf" / f"{h}.white" for h in ("lh", "rh")]), adult_tal)
    tree = cKDTree(adult_mni[::10])
    tal = {}
    for child in cfg["anatomy"]["school"]:
        name = child["subject"]
        scanner = mni_cortex([root / name / "surf" / f"{h}.white" for h in ("lh", "rh")])
        for dest, who in tal_files.items():
            if not (dest.startswith(f"{name}/") or who == name):
                continue
            c = to_mni(scanner, fsio.read_talairach(root / dest))
            tal[f"{child['key']}/{dest}"] = dict(
                file_owner=who, own_file=who == name, stored_in_own_folder=dest.startswith(f"{name}/"),
                centroid_offset_from_adult_mm=float(np.linalg.norm(c.mean(axis=0) - adult_mni.mean(axis=0))),
                median_distance_to_adult_cortex_mm=float(np.median(tree.query(c[::10])[0])))
    out["talairach_files"] = tal
    io.write_json(out, OUT / "school_anatomy_checks.json")
    print(json.dumps({k: v for k, v in out.items() if k not in ("status", "upper_head_definition")}, indent=1))


if __name__ == "__main__":
    main()
