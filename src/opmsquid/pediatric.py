"""Pediatric fixed-helmet versus head-adaptive comparison (G3B): head size, source-blind placements
of a head in the fixed adult Neuromag helmet, a counterfactual reduced-mismatch helmet, and
summaries over homologous sources.

Frames: the Neuromag device frame (+x right, +y anterior, +z up) and each subject's head frame
(Neuromag convention from its own fiducials). A rigid head motion ``m`` (device frame) changes the
device-to-head transform to ``dev_head @ inv(m)`` (as in ``g2.head_position_variants``).

Placement rules (A-G3-PLACE; chosen without reference to any test source):
* ``centred``: the head frame where the adult's is in the sample recording (the adult's measured
  device-to-head transform): the ear line stays where the adult's was, and a smaller head has a
  wider gap at the vertex.
* contact: from a start pose, the head is moved along one device direction until the nearest
  magnetometer coil centre is ``clearance`` (20 mm = the 18-mm Dewar spacing HW-18mm + 2 mm) from
  the scalp; ``top`` moves up (+z, the usual pediatric positioning against the top of the helmet),
  ``back`` moves posterior (-y, the occiput against the back of the helmet, as in a supine child).
* bounded variants (``x+5mm`` ... ``roll-5deg``): the centred head translated by +-5 mm along device
  x and y, pitched by +-10 deg or rolled by +-5 deg about the head origin, then raised to the same
  top contact where there is room (a pose already within the clearance, e.g. an adult head shifted
  towards the helmet wall, stays as it is). For a child these are variants of ``top``.
* ``x-centred``: the centred head shifted along device x until the median magnetometer-to-scalp
  distances of the left and right helmet halves are equal (a lateral centring that uses geometry
  only; the adult's measured pose leaves a head with other fiducials off-centre), then top contact.
* ``top-18mm``: the centred head raised until the nearest coil is at the 18-mm Dewar spacing itself
  (true contact).
Feasibility: every magnetometer coil centre at least ``dewar`` (18 mm) from the scalp.
"""
from __future__ import annotations

import numpy as np
import mne
from scipy.spatial import ConvexHull, cKDTree

from . import neuromag

CLEARANCE = 0.020  # m, contact rule: nearest magnetometer coil centre to scalp (HW-18mm + 2 mm), as G2's 'well_fitted'
DEWAR = 0.018  # m, feasibility: no magnetometer coil centre closer to the scalp (HW-18mm)


# ----------------------------------------------------------------------------------------------
# head size
def scalp_head_frame(subject) -> np.ndarray:
    mh = np.linalg.inv(subject.trans["trans"])
    return subject.scalp.rr @ mh[:3, :3].T + mh[:3, 3]


def head_size(subject, step: float = 0.002, half_width: float = 0.0015) -> dict:
    """Head dimensions [mm] from the dense scalp, head frame (z = 0: the plane of LPA, RPA and
    nasion): occipitofrontal circumference (largest perimeter of the convex hull of scalp sections
    parallel to that plane, 0-60 mm above it), maximum breadth and length and the vertex height
    above the plane, the convex-hull volume of the scalp above the plane [cm^3] and the
    inter-auricular (LPA-RPA) distance."""
    rr = scalp_head_frame(subject)
    cap = rr[rr[:, 2] >= 0]
    best, z_best = 0.0, 0.0
    for z in np.arange(0.0, 0.060, step):
        sec = rr[np.abs(rr[:, 2] - z) < half_width][:, :2]
        if len(sec) < 10:
            continue
        hull = ConvexHull(sec)
        if hull.area > best:  # in 2-D, ConvexHull.area is the perimeter
            best, z_best = hull.area, z
    fids = subject.fiducials
    if fids is None:
        from .anatomy import _sample_fiducials

        fids = _sample_fiducials()
    return dict(ofc_mm=best * 1e3, ofc_plane_height_mm=z_best * 1e3, breadth_mm=float(np.ptp(cap[:, 0]) * 1e3),
                length_mm=float(np.ptp(cap[:, 1]) * 1e3), vertex_height_mm=float(cap[:, 2].max() * 1e3),
                cap_volume_cm3=float(ConvexHull(cap).volume * 1e6),
                inter_auricular_mm=float(np.linalg.norm(fids["lpa"] - fids["rpa"]) * 1e3))


# ----------------------------------------------------------------------------------------------
# helmet geometry and placements
def translate(d) -> np.ndarray:
    m = np.eye(4)
    m[:3, 3] = d
    return m


def rotate_about(axis: int, deg: float, point) -> np.ndarray:
    """Rotation by ``deg`` about device axis ``axis`` (0 x, 1 y, 2 z) through ``point`` (device frame)."""
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    i, j = [k for k in range(3) if k != axis]
    r = np.eye(4)
    r[i, i], r[i, j], r[j, i], r[j, j] = c, -s, s, c
    return translate(point) @ r @ translate(-np.asarray(point))


def moved(dev_head: np.ndarray, motion: np.ndarray) -> np.ndarray:
    return dev_head @ np.linalg.inv(motion)


def magnetometer_positions(info: mne.Info) -> np.ndarray:
    """Magnetometer coil centres, device frame."""
    kinds = neuromag.channel_kinds(info)
    return np.array([c["loc"][:3] for c, k in zip(info["chs"], kinds) if k == "mag"])


class HelmetFit:
    """Distances between the Neuromag magnetometer coils of ``info`` and a subject's scalp."""

    def __init__(self, info: mne.Info, subject):
        from mne.surface import _CheckInside

        self.coils = magnetometer_positions(info)
        self.tree = cKDTree(scalp_head_frame(subject))
        skin = next(s for s in subject.bem_surfaces if s["id"] == mne.io.constants.FIFF.FIFFV_BEM_SURF_ID_HEAD)
        self._inside_mri = _CheckInside(skin)
        self._head_mri = subject.trans["trans"]

    def head_coils(self, dev_head):
        return self.coils @ dev_head[:3, :3].T + dev_head[:3, 3]

    def distances(self, dev_head) -> np.ndarray:
        return self.tree.query(self.head_coils(dev_head))[0]

    def any_inside(self, dev_head) -> bool:
        p = self.head_coils(dev_head)
        p_mri = p @ self._head_mri[:3, :3].T + self._head_mri[:3, 3]
        return bool(np.any(self._inside_mri(p_mri)))

    def contact(self, dev_head: np.ndarray, direction, clearance: float = CLEARANCE, step: float = 0.0005,
                max_move: float = 0.100) -> tuple[np.ndarray, float]:
        """Move the head along the device-frame unit ``direction`` from ``dev_head`` until the next
        step would bring the nearest coil closer than ``clearance`` to the scalp. A pose already
        closer than that (e.g. after a sideways shift towards the helmet wall) is not moved; its
        feasibility is judged by the Dewar spacing (``describe``). Returns the transform and the
        distance moved [m]."""
        u = np.asarray(direction, float) / np.linalg.norm(direction)
        if self.distances(dev_head).min() <= clearance:
            return dev_head, 0.0
        prev = 0.0
        for k in range(1, int(max_move / step) + 1):
            d = k * step
            if self.distances(moved(dev_head, translate(d * u))).min() <= clearance:
                return moved(dev_head, translate(prev * u)), prev
            prev = d
        raise RuntimeError("no contact within the search range")

    def describe(self, dev_head) -> dict:
        d = self.distances(dev_head)
        return dict(min_dist_mm=float(d.min() * 1e3), median_dist_mm=float(np.median(d) * 1e3), max_dist_mm=float(d.max() * 1e3),
                    feasible=bool(d.min() >= DEWAR and not self.any_inside(dev_head)))


def lateral_centring(info: mne.Info, subject, dev_head: np.ndarray, max_shift: float = 0.020,
                     step: float = 0.0005) -> tuple[np.ndarray, float]:
    """``dev_head`` with the head shifted along device x (within +-``max_shift``) so that the median
    magnetometer-to-scalp distance of the left helmet half (MNE's left selections) equals that of the
    right half. Returns the transform and the shift [m] (positive: towards device +x, the right)."""
    fit = HelmetFit(info, subject)
    regions = helmet_regions(info)
    left = np.concatenate([v for k, v in regions.items() if k.startswith("left")])
    right = np.concatenate([v for k, v in regions.items() if k.startswith("right")])
    best = None
    for dx in np.arange(-max_shift, max_shift + step / 2, step):
        t = moved(dev_head, translate([dx, 0.0, 0.0]))
        d = fit.distances(t)
        imbalance = abs(float(np.median(d[left]) - np.median(d[right])))
        if best is None or imbalance < best[0] - 1e-12:
            best = (imbalance, float(dx), t)
    return best[2], best[1]


def placements(info: mne.Info, subject, base_dev_head: np.ndarray, translation: float = 0.005, pitch_deg: float = 10.0,
               roll_deg: float = 5.0, clearance: float = CLEARANCE) -> dict:
    """Source-blind head placements in the helmet of ``info`` (see the module docstring): name ->
    dict(trans=device-to-head 4x4, rule, distance stats, feasible)."""
    fit = HelmetFit(info, subject)
    origin = np.linalg.inv(base_dev_head)[:3, 3]  # head origin, device frame
    out = {"centred": dict(trans=base_dev_head, rule="adult measured device-to-head transform", moved_mm=0.0)}
    t, d = fit.contact(base_dev_head, (0, 0, 1), clearance)
    out["top"] = dict(trans=t, rule=f"centred, then up (+z) to {clearance * 1e3:g}-mm contact", moved_mm=d * 1e3)
    t, d = fit.contact(base_dev_head, (0, -1, 0), clearance)
    out["back"] = dict(trans=t, rule=f"centred, then posterior (-y) to {clearance * 1e3:g}-mm contact", moved_mm=d * 1e3)
    variants = {f"x{s}{translation * 1e3:g}mm": translate(sg * translation * np.eye(3)[0]) for s, sg in (("+", 1), ("-", -1))}
    variants.update({f"y{s}{translation * 1e3:g}mm": translate(sg * translation * np.eye(3)[1]) for s, sg in (("+", 1), ("-", -1))})
    variants.update({f"pitch{s}{pitch_deg:g}deg": rotate_about(0, sg * pitch_deg, origin) for s, sg in (("+", 1), ("-", -1))})
    variants.update({f"roll{s}{roll_deg:g}deg": rotate_about(1, sg * roll_deg, origin) for s, sg in (("+", 1), ("-", -1))})
    for name, m in variants.items():
        t, d = fit.contact(moved(base_dev_head, m), (0, 0, 1), clearance)
        out[name] = dict(trans=t, rule=f"centred, {name}, then up to contact where there is room", moved_mm=d * 1e3)
    pose, dx = lateral_centring(info, subject, base_dev_head)
    t, d = fit.contact(pose, (0, 0, 1), clearance)
    out["x-centred"] = dict(trans=t, pose=pose, shift_x_mm=dx * 1e3, moved_mm=d * 1e3,
                            rule=f"centred, shifted {dx * 1e3:+.1f} mm along device x to equal left/right median gaps, then up to "
                                 f"{clearance * 1e3:g}-mm contact")
    t, d = fit.contact(base_dev_head, (0, 0, 1), DEWAR)
    out["top-18mm"] = dict(trans=t, moved_mm=d * 1e3, rule=f"centred, then up to {DEWAR * 1e3:g}-mm contact (true contact)")
    for v in out.values():
        v.update(fit.describe(v["trans"]))
    return out


def with_dev_head(info: mne.Info, dev_head: np.ndarray) -> mne.Info:
    out = info.copy()
    with out._unlock():
        out["dev_head_t"] = mne.transforms.Transform("meg", "head", np.asarray(dev_head, float))
    return out


def scaled_helmet(info: mne.Info, k: float, centre) -> mne.Info:
    """Counterfactual helmet: every coil centre scaled by ``k`` about ``centre`` (device frame);
    coil orientations, sizes, integration rules and noise unchanged."""
    c = np.asarray(centre, float)
    out = info.copy()
    with out._unlock():
        for ch in out["chs"]:
            ch["loc"][:3] = c + k * (np.asarray(ch["loc"][:3]) - c)
    return out


def counterfactual_helmet(info: mne.Info, subject, base_dev_head: np.ndarray, k: float, step: float = 0.005) -> dict:
    """Reduced head-helmet mismatch (A-G3-COUNTERFACTUAL; a mechanistic control, not a pediatric
    SQUID system): the helmet scaled by ``k`` (the child/adult head-circumference ratio) about the
    head origin of the ``centred`` placement, which maps the adult's measured fit onto a
    geometrically similar head exactly (every gap scaled by k). If that brings a coil closer to the
    scalp than the Dewar spacing, k is increased in ``step`` until it does not."""
    centre = np.linalg.inv(base_dev_head)[:3, 3]
    k0 = float(k)
    while True:
        inf = scaled_helmet(info, k, centre)
        fit = HelmetFit(inf, subject)
        desc = fit.describe(base_dev_head)
        if desc["feasible"] or k >= 1.0:
            break
        k += step
    return dict(k=float(k), k_nominal=k0, info=inf, trans=base_dev_head,
                rule=f"helmet scaled by {k:.3f} about the head origin (head-circumference ratio {k0:.3f}); centred placement",
                **desc)


def helmet_regions(info: mne.Info) -> dict:
    """MNE's Vectorview channel selections (Left/Right temporal, parietal, occipital, frontal) as
    magnetometer index sets (order of ``magnetometer_positions``)."""
    kinds = neuromag.channel_kinds(info)
    mag_names = [c["ch_name"] for c, k in zip(info["chs"], kinds) if k == "mag"]
    out = {}
    for side in ("Left", "Right"):
        for region in ("frontal", "temporal", "parietal", "occipital"):
            names = set(mne.read_vectorview_selection([f"{side}-{region}"], info=info))
            out[f"{side.lower()} {region}"] = np.array([i for i, n in enumerate(mag_names) if n in names])
    return out


# ----------------------------------------------------------------------------------------------
# summaries
def weighted_median(x: np.ndarray, w: np.ndarray | None = None) -> float:
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    if w is None:
        return float(np.median(x[ok])) if ok.any() else float("nan")
    w = np.asarray(w, float)[ok]
    x = x[ok]
    if not len(x):
        return float("nan")
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return float(x[o][np.searchsorted(c, 0.5 * c[-1])])


def grouped_bootstrap(x: np.ndarray, groups: np.ndarray, w: np.ndarray | None, rng: np.random.Generator,
                      n_boot: int = 1000) -> list[float]:
    """95 % interval of the (area-weighted) median of ``x`` when whole groups (parcels) are
    resampled with replacement."""
    labels = np.unique(groups)
    idx = [np.flatnonzero(groups == g) for g in labels]
    w = np.ones(len(x)) if w is None else np.asarray(w, float)
    boot = []
    for _ in range(n_boot):
        take = np.concatenate([idx[j] for j in rng.integers(0, len(idx), len(idx))])
        boot.append(weighted_median(x[take], w[take]))
    return [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]


def target_areas(cortex, targets: np.ndarray) -> np.ndarray:
    """Cortical area [m^2] represented by each target: usable full-resolution vertex areas summed
    over each target's nearest-target (Euclidean) cell (area weighting of summaries)."""
    usable = np.flatnonzero(cortex.usable)
    owner = cKDTree(cortex.rr[targets]).query(cortex.rr[usable])[1]
    return np.bincount(owner, weights=cortex.area[usable], minlength=len(targets))
