"""Real Neuromag sensor geometry: the SQUID comparators of G1B, G1C, G2 and G3.

Geometry source: the MNE sample recording ``sample_audvis_raw.fif`` (MGH Elekta Neuromag
Vectorview, 306 channels = 102 sensor locations, each with one magnetometer and two orthogonal
planar gradiometers), with its measured device-to-head transform and the ``sample`` subject's
head-to-MRI transform. This is representative Vectorview/TRIUX helmet geometry, not the MRN
installation (which has no geometry file available here).

Coil types: MRN's TRIUX uses T3 sensors, planar gradiometer 3014 and magnetometer 3024
(www.mrn.org/collaborate/elekta-neuromag-meg, read 2026-09-30). The sample file stores 3012
(planar gradiometer T1) and 3024. In MNE 1.13.2's ``coil_def.dat``, 3012 and 3014 share size
(26.39 mm), baseline (16.80 mm) and the normal/accurate integration points; only the 1-point
representation differs (0.3 mm offset). ``load_info(variant="T3")`` therefore sets 3014 on all
gradiometers; ``variant="file"`` keeps the stored types. Magnetometers are never relabelled
as gradiometers or vice versa.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import mne
from mne.io.constants import FIFF

from . import paths

T3_COILS = {"grad": int(FIFF.FIFFV_COIL_VV_PLANAR_T3), "mag": int(FIFF.FIFFV_COIL_VV_MAG_T3)}  # 3014, 3024
GRAD_TYPES = {int(FIFF.FIFFV_COIL_VV_PLANAR_T1), int(FIFF.FIFFV_COIL_VV_PLANAR_T2), int(FIFF.FIFFV_COIL_VV_PLANAR_T3)}
MAG_TYPES = {int(FIFF.FIFFV_COIL_VV_MAG_T1), int(FIFF.FIFFV_COIL_VV_MAG_T2), int(FIFF.FIFFV_COIL_VV_MAG_T3)}
RAW_FILE = "sample_audvis_raw.fif"
TRANS_FILE = "sample_audvis_raw-trans.fif"


def load_info(variant: str = "T3") -> mne.Info:
    """MEG-only measurement info of the sample recording: 306 channels, no bad-channel
    exclusion (geometry only), no projectors or compensation.

    variant: 'T3' (3014/3024, MNE 'accurate' integration), 'file' (stored types) or 'T3-4pt'
    (study coils 9014/9024 = T3 with MNE's 4-point 'normal' rule, as in Hunold et al. 2016;
    forward computations then need ``mne.use_coil_def(opm.coil_def_file())``)."""
    if variant == "T3-4pt":
        from .opm import FOUR_POINT_COILS

        info = load_info("T3")
        with info._unlock():
            for ch in info["chs"]:
                ch["coil_type"] = FOUR_POINT_COILS[int(ch["coil_type"])]
        return info
    raw_info = mne.io.read_info(paths.require(paths.SAMPLE_MEG / RAW_FILE, "MNE sample recording"), verbose=False)
    info = mne.pick_info(raw_info, mne.pick_types(raw_info, meg=True, eeg=False, stim=False, eog=False, exclude=[]))
    with info._unlock():
        info["bads"] = []
        info["projs"] = []
        info["comps"] = []
        for ch in info["chs"]:
            ct = int(ch["coil_type"])
            if ct not in GRAD_TYPES | MAG_TYPES:
                raise ValueError(f"unexpected MEG coil type {ct} for {ch['ch_name']}")
            if variant == "T3":
                ch["coil_type"] = T3_COILS["grad"] if ct in GRAD_TYPES else T3_COILS["mag"]
            elif variant != "file":
                raise ValueError("variant must be 'T3' or 'file'")
    return info


def channel_kinds(info: mne.Info) -> np.ndarray:
    """'grad' (unit T/m) or 'mag' (unit T) for every MEG channel, independent of coil type."""
    return np.array(["grad" if ch["unit"] == FIFF.FIFF_UNIT_T_M else "mag" for ch in info["chs"]
                     if ch["kind"] == FIFF.FIFFV_MEG_CH])


def head_mri_trans() -> mne.transforms.Transform:
    return mne.read_trans(paths.require(paths.SAMPLE_MEG / TRANS_FILE, "sample head-MRI transform"))


@dataclass
class SensorGeometry:
    names: list
    kind: np.ndarray  # "mag" or "grad"
    pos: np.ndarray  # (n, 3) coil centres [m] in `frame`
    normal: np.ndarray  # (n, 3) coil normals (unit)
    site: np.ndarray  # (n,) index of the 102 sensor locations (triplets)
    frame: str


def sensor_geometry(info: mne.Info, frame: str = "head", trans: mne.transforms.Transform | None = None
                    ) -> SensorGeometry:
    """Coil centres and normals of all MEG channels in the head or MRI frame, with the index
    of the sensor location each channel belongs to (channels whose centres coincide)."""
    chs = [ch for ch in info["chs"] if ch["kind"] == FIFF.FIFFV_MEG_CH]
    pos = np.array([ch["loc"][:3] for ch in chs])
    nrm = np.array([ch["loc"][9:12] for ch in chs])
    t = info["dev_head_t"]["trans"]
    if frame == "mri":
        t = (trans if trans is not None else head_mri_trans())["trans"] @ t
    elif frame != "head":
        raise ValueError("frame must be 'head' or 'mri'")
    pos = pos @ t[:3, :3].T + t[:3, 3]
    nrm = nrm @ t[:3, :3].T
    kind = channel_kinds(info)
    # sensor locations: channels of one triplet share the coil centre (to < 0.1 mm)
    mags = np.flatnonzero(kind == "mag")
    site = np.argmin(np.linalg.norm(pos[:, None, :] - pos[None, mags, :], axis=2), axis=1)
    return SensorGeometry([ch["ch_name"] for ch in chs], kind, pos, nrm, site, frame)


def ray_mesh_distance(origins: np.ndarray, directions: np.ndarray, rr: np.ndarray, tris: np.ndarray,
                      chunk: int = 64) -> np.ndarray:
    """Distance along each ray (origin + t direction, t > 0) to its first intersection with a
    triangle mesh (Moller-Trumbore); NaN if the ray misses."""
    v0, v1, v2 = rr[tris[:, 0]], rr[tris[:, 1]], rr[tris[:, 2]]
    e1, e2 = v1 - v0, v2 - v0
    out = np.full(len(origins), np.nan)
    for start in range(0, len(origins), chunk):
        o = origins[start:start + chunk, None, :]
        d = directions[start:start + chunk, None, :] / np.linalg.norm(directions[start:start + chunk], axis=1)[:, None, None]
        p = np.cross(d, e2[None])
        det = np.sum(e1[None] * p, axis=2)
        ok = np.abs(det) > 1e-14
        inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
        tvec = o - v0[None]
        u = np.sum(tvec * p, axis=2) * inv
        q = np.cross(tvec, e1[None])
        v = np.sum(d * q, axis=2) * inv
        t = np.sum(e2[None] * q, axis=2) * inv
        hit = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 0)
        t = np.where(hit, t, np.inf)
        best = t.min(axis=1)
        out[start:start + chunk] = np.where(np.isfinite(best), best, np.nan)
    return out
