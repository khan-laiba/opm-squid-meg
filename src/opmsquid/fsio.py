"""FreeSurfer files without nibabel: triangle surfaces (with their volume information) and the
``talairach.xfm`` transform, and MNE's surface source space for surfaces read this way.

The study's environment has no nibabel, which MNE needs to read FreeSurfer geometry. These readers
follow the FreeSurfer file formats as ``nibabel.freesurfer.io`` reads them; ``surface_source_space``
repeats the steps of MNE 1.13's ``setup_source_space`` on surfaces already in memory. Both are
checked against MNE's stored source space of the sample subject (``tests/test_fsio.py``). Used to
prepare the school-aged children (``scripts/prepare_school_subjects.py``).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

TRIANGLE_MAGIC = 16777214
_VOLUME_KEYS = ("valid", "filename", "volume", "voxelsize", "xras", "yras", "zras", "cras")


def _fread3(fobj) -> int:
    b = np.fromfile(fobj, ">u1", 3).astype(np.int64)
    if len(b) != 3:
        raise ValueError("truncated FreeSurfer file")
    return int((b[0] << 16) + (b[1] << 8) + b[2])


def read_geometry(fname, read_metadata: bool = False):
    """Vertices [mm] (float64) and triangles (int64) of a FreeSurfer triangle surface; with
    ``read_metadata`` also its volume information (dict, empty if the file has none)."""
    with open(fname, "rb") as f:
        magic = _fread3(f)
        if magic != TRIANGLE_MAGIC:
            raise ValueError(f"{fname}: not a FreeSurfer triangle surface (magic {magic})")
        f.readline()  # creation stamp
        f.readline()
        vnum, fnum = (int(x) for x in np.fromfile(f, ">i4", 2))
        rr = np.fromfile(f, ">f4", vnum * 3)
        tris = np.fromfile(f, ">i4", fnum * 3)
        if rr.size != vnum * 3 or tris.size != fnum * 3:
            raise ValueError(f"{fname}: truncated surface")
        rr = rr.reshape(vnum, 3).astype(np.float64)
        tris = tris.reshape(fnum, 3).astype(np.int64)
        info = _read_volume_info(f) if read_metadata else None
    if tris.min() < 0 or tris.max() >= vnum:
        raise ValueError(f"{fname}: triangle index out of range")
    return (rr, tris, info) if read_metadata else (rr, tris)


def _read_volume_info(f) -> dict:
    head = np.fromfile(f, ">i4", 1)
    if len(head) == 0:
        return {}
    if not np.array_equal(head, [20]):
        head = np.concatenate([head, np.fromfile(f, ">i4", 2)])
        if not (np.array_equal(head, [2, 0, 20]) or np.array_equal(head, [2, 1, 20])):
            return {}
    out = {}
    for key in _VOLUME_KEYS:
        name, _, value = f.readline().decode("utf-8").partition("=")
        if name.strip() != key:
            raise ValueError(f"volume information: expected '{key}', got '{name.strip()}'")
        value = value.strip()
        out[key] = value if key in ("valid", "filename") else np.array(value.split(), int if key == "volume" else float)
    return out


def surface_ras_to_scanner(info: dict) -> np.ndarray:
    """4 x 4 transform [mm] from FreeSurfer surface RAS (MNE's MRI frame) to scanner RAS, from the
    volume information of a surface file: Norig inv(Torig) of the volume the surface was made from
    (as MNE computes it from the MRI header)."""
    dims, vs = np.asarray(info["volume"], float), np.asarray(info["voxelsize"], float)
    m = np.column_stack([info["xras"], info["yras"], info["zras"]]).astype(float) * vs
    norig = np.eye(4)
    norig[:3, :3] = m
    norig[:3, 3] = np.asarray(info["cras"], float) - m @ (dims / 2.0)
    torig = np.array([[-vs[0], 0.0, 0.0, vs[0] * dims[0] / 2.0], [0.0, 0.0, vs[2], -vs[2] * dims[2] / 2.0],
                      [0.0, -vs[1], 0.0, vs[1] * dims[1] / 2.0], [0.0, 0.0, 0.0, 1.0]])
    return norig @ np.linalg.inv(torig)


def read_talairach(fname) -> np.ndarray:
    """4 x 4 linear transform [mm] of a FreeSurfer ``talairach.xfm`` (scanner RAS -> MNI305)."""
    lines = Path(fname).read_text().splitlines()
    i = next((k for k, line in enumerate(lines) if line.startswith("Linear_Transform")), None)
    if i is None:
        raise ValueError(f"{fname}: no Linear_Transform")
    rows = [[float(x) for x in line.strip().rstrip(";").split()] for line in lines[i + 1:i + 4]]
    if len(rows) != 3 or any(len(r) != 4 for r in rows):
        raise ValueError(f"{fname}: malformed Linear_Transform")
    return np.vstack([np.array(rows), [0.0, 0.0, 0.0, 1.0]])


def mni_fiducials(surface_info: dict, talairach: np.ndarray) -> dict:
    """LPA, nasion and RPA in the subject's MRI frame [m]: fsaverage's fiducials (MNE's package
    data, MNI305 [m]) mapped through the inverse of the subject's MRI -> MNI305 transform, as
    ``mne.coreg.get_mni_fiducials`` does."""
    import mne
    from mne.io.constants import FIFF

    fids, frame = mne.io.read_fiducials(Path(mne.__file__).parent / "data" / "fsaverage" / "fsaverage-fiducials.fif",
                                        verbose=False)
    if frame != FIFF.FIFFV_COORD_MRI:
        raise ValueError("fsaverage fiducials are not in the MRI frame")
    mni_mri = np.linalg.inv(talairach @ surface_ras_to_scanner(surface_info))  # mm
    key = {FIFF.FIFFV_POINT_LPA: "lpa", FIFF.FIFFV_POINT_NASION: "nasion", FIFF.FIFFV_POINT_RPA: "rpa"}
    return {key[f["ident"]]: (mni_mri @ np.r_[np.asarray(f["r"], float) * 1e3, 1.0])[:3] / 1e3
            for f in fids if f["ident"] in key}


def surface_source_space(whites, spheres, subject: str, spacing: str = "oct6"):
    """MNE's surface source space (``setup_source_space(subject, spacing, surface='white',
    add_dist='patch')``) for surfaces already read: ``whites`` and ``spheres`` hold (rr [mm], tris)
    for lh and rh. The same steps as MNE 1.13 (the sphere's nearest vertices to the subdivided
    octahedron or icosahedron, displaced to a free neighbour on double occupation; then the fill-in
    and the patch statistics), so the vertices are MNE's."""
    from mne.io.constants import FIFF
    from mne.source_space import _source_space as S

    stype, _, ico_surf, _ = S._check_spacing(spacing)
    if stype not in ("ico", "oct"):
        raise ValueError("only ico/oct spacings are supported")
    ico_rr = np.array(ico_surf["rr"], float)
    S._normalize_vectors(ico_rr)
    out = []
    for (rr, tris), (srr, _), sid in zip(whites, spheres, (FIFF.FIFFV_MNE_SURF_LEFT_HEMI, FIFF.FIFFV_MNE_SURF_RIGHT_HEMI)):
        surf = dict(rr=np.array(rr, float), tris=np.array(tris), ntri=len(tris), use_tris=np.array(tris), np=len(rr))
        S.complete_surface_info(surf, False, copy=False)
        sphere = np.array(srr, float)
        S._normalize_vectors(sphere)
        if len(sphere) != surf["np"]:
            raise RuntimeError("white and sphere surfaces differ in their number of vertices")
        mmap = S._compute_nearest(sphere, ico_rr)
        surf["inuse"] = np.zeros(surf["np"], int)
        for k in range(len(mmap)):
            if surf["inuse"][mmap[k]]:
                neigh = S._get_surf_neighbors(surf, mmap[k])
                free = np.where(np.logical_not(surf["inuse"][neigh]))[0]
                if len(free) == 0:
                    raise RuntimeError(f"no free neighbour for vertex {k} / {len(mmap)}")
                mmap[k] = neigh[free[-1]]
            elif mmap[k] < 0 or mmap[k] > surf["np"]:
                raise RuntimeError(f"map number out of range ({mmap[k]})")
            surf["inuse"][mmap[k]] = True
        surf["use_tris"] = np.array([mmap[t] for t in ico_surf["tris"]], np.int32)
        surf["nuse_tri"] = len(surf["use_tris"])
        surf["nuse"] = np.sum(surf["inuse"])
        surf["vertno"] = np.where(surf["inuse"])[0]
        sizes = S._normalize_vectors(surf["nn"])
        surf["inuse"][sizes <= 0] = False
        surf["nuse"] = np.sum(surf["inuse"])
        surf["subject_his_id"] = subject
        surf.update(dist=None, dist_limit=None, nearest=None, type="surf", nearest_dist=None, pinfo=None, patch_inds=None, id=sid,
                    coord_frame=FIFF.FIFFV_COORD_MRI)
        surf["rr"] /= 1000.0
        for k in ("tri_area", "tri_cent", "tri_nn", "neighbor_tri"):
            surf.pop(k, None)
        out.append(surf)
    src = S.SourceSpaces(out, dict(working_dir=os.getcwd(), command_line=f"opmsquid.fsio.surface_source_space({subject}, {spacing})"))
    S.add_source_space_distances(src, dist_limit=0.0, verbose=False)
    return src
