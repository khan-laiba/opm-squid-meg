"""Realistic adult OPM-Neuromag comparison (G2): arrays, sources, noise models, metrics.

Arrays (sample subject, measured head position in the Neuromag helmet):
  squid        Neuromag T3, 102 magnetometers + 204 planar gradiometers (comparators: 'mag',
               'grad', 'combined')
  opm_matched  matched-site OPM array (coverage control; the Neuromag sites that fit the OPM
               placement rules, 95 of 102 on the sample head with the v3 clearance rule)
  opm204       204 sites spread evenly over the dense array (channel-budget control vs the 204
               gradiometers)
  opm_dense    a dense single-axis OPM array under the 17-mm packing rule, greedy farthest-point
               construction ("full system"; 205 sites on the sample head with the v3 clearance rule). Not proven maximal;
               306 single-axis channels appear infeasible on this head (A-OPM-PACK).
OPM sensitive axes follow the smooth BEM head-surface normal (A-OPM-AXIS). Sites avoid the ears,
the ear pinna and the edge of the MRI field of view (A-OPM-COVER); a package may sit at most 5 mm
beyond its nominal standoff to clear the scalp (A-OPM-CLEAR).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import mne
from mne.io.constants import FIFF

from . import anatomy, background, environment, forward, goldenholz, metrics, neuromag, noisemodel, opm, paths

SQUID_ASD = {"mag": 3.5e-15, "grad": 3.6e-13}  # HW-mag-noise, HW-grad-noise (typical white noise)
OPM_ASD_SWEEP = (7e-15, 10e-15, 15e-15, 20e-15, 30e-15)  # A-OPM-NOISE
BEM_CONDUCTIVITY = (0.3, 0.006, 0.3)  # MNE default 3-layer (scalp, skull, brain)


@dataclass
class Array:
    name: str
    info: mne.Info
    kinds: np.ndarray  # 'mag' | 'grad' per channel
    coil_def: Path | None
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.kinds)


def build_arrays(subject: anatomy.Subject, digitisation: mne.Info, scalp_gap: float = 0.0,
                 dev_head_t: np.ndarray | None = None) -> dict[str, Array]:
    squid_info = neuromag.load_info("T3")
    if dev_head_t is not None:
        with squid_info._unlock():
            squid_info["dev_head_t"] = mne.transforms.Transform("meg", "head", dev_head_t)
    arrays = {"squid": Array("squid", squid_info, neuromag.channel_kinds(squid_info), None,
                             dict(sites=102, channels=306, axes="1 mag + 2 planar grad per site"))}
    arrays["opm_matched"] = matched_opm(subject, digitisation, scalp_gap)
    for name in DENSE:
        arrays[name] = dense_opm(subject, digitisation, name, scalp_gap)
    return arrays


# dense single-axis arrays: (farthest-point scalp spacing [m] of the densest feasible array,
# number of sites taken from it by farthest-point sampling (None = all), role)
DENSE = {"opm204": (0.015, 204, "channel-budget control vs 204 gradiometers"),
         "opm_dense": (0.015, None, "densest feasible single-axis array (full system)")}


AXIS_RADIUS = 0.015  # A-OPM-AXIS: sensitive axis = BEM head-surface normal averaged within 15 mm
EAR_CLEARANCE = 0.020  # A-OPM-COVER: no site within 20 mm of the preauricular points


def skin_surface(subject: anatomy.Subject) -> anatomy.Surface:
    """The smooth BEM head surface (MRI frame, outward normals), used for OPM sensitive axes."""
    return anatomy._outward(next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD))


def matched_opm(subject: anatomy.Subject, digitisation: mne.Info, scalp_gap: float = 0.0,
                squid_info: mne.Info | None = None) -> Array:
    """Matched-site array: the Neuromag sites of ``squid_info`` (default: the sample recording, i.e.
    the measured head position) projected onto the scalp, under the OPM placement rules."""
    skin = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    squid_info = neuromag.load_info("T3") if squid_info is None else squid_info
    arr, rep = opm.matched_to_neuromag(squid_info, subject.trans, subject.scalp, digitisation,
                                       scalp_gap=scalp_gap, normal_radius=AXIS_RADIUS, axis_surface=skin_surface(subject),
                                       ear_clearance=EAR_CLEARANCE, outer_skin=skin)
    return _opm_array("opm_matched", arr, dict(role="matched-site coverage control", **_rep(rep)))


def dense_opm(subject: anatomy.Subject, digitisation: mne.Info, name: str, scalp_gap: float = 0.0) -> Array:
    """Densest feasible single-axis array (17-mm packing rule), or an evenly spread subset of it
    with a fixed channel count (farthest-point sampling of the sensing centres)."""
    spacing, n_sites, role = DENSE[name]
    skin = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    arr, rep = opm.dense_array(subject.scalp, subject.trans, digitisation, spacing, scalp_gap=scalp_gap, outer_skin=skin,
                               normal_radius=AXIS_RADIUS, axis_surface=skin_surface(subject), ear_clearance=EAR_CLEARANCE)
    if n_sites is not None:
        if n_sites > len(arr.pos):
            raise ValueError(f"{name}: only {len(arr.pos)} feasible sites")
        arr = opm.subset(arr, opm.farthest_point_subset(arr.pos, n_sites), f"dense OPM subset ({n_sites} sites)")
        parent = {f"parent_{k}": rep[k] for k in ("n_moved_out", "max_extra_shift_mm") if k in rep}
        rep = {k: v for k, v in rep.items() if k not in ("n_moved_out", "max_extra_shift_mm")}
        rep = dict(rep, n_sites=n_sites, min_spacing_mm=float(opm.min_spacing(arr.pos).min() * 1e3),
                   median_spacing_mm=float(np.median(opm.min_spacing(arr.pos)) * 1e3), subset_of=rep["n_sites"], **parent)
    return _opm_array(name, arr, dict(role=role, **rep))


def _rep(rep):
    return dict(n_sites=rep["n_kept"], excluded_sites=rep["excluded_sites"], extra_shift_mm=rep["extra_shift_mm"],
                min_spacing_mm=float(np.min(rep["min_spacing_mm"])))


def _opm_array(name, arr, meta) -> Array:
    info = opm.make_info(arr)
    meta.update(channels=len(arr.pos), axes="1 (scalp normal)", standoff_mm=arr.standoff * 1e3,
                scalp_gap_mm=arr.scalp_gap * 1e3)
    return Array(name, info, np.array(["mag"] * len(arr.pos)), opm.coil_def_file(), meta)


def with_scalp_gap(array: Array, gap: float) -> Array:
    """The same OPM sites moved outward along their sensitive axes by ``gap`` [m] (A-OPM-GAP):
    the helmet sits farther from the scalp and every sensor keeps its site, so a gap variant
    differs from the primary array only by the gap (clearance can only increase)."""
    info = array.info.copy()
    with info._unlock():
        for ch in info["chs"]:
            ch["loc"][:3] = ch["loc"][:3] + gap * ch["loc"][9:12]
    meta = dict(array.meta, scalp_gap_mm=array.meta.get("scalp_gap_mm", 0.0) + gap * 1e3)
    return Array(array.name, info, array.kinds.copy(), array.coil_def, meta)


def channel_sets(array: Array) -> dict[str, np.ndarray]:
    if array.name == "squid":
        return {"mag": array.kinds == "mag", "grad": array.kinds == "grad", "combined": np.ones(array.n, bool)}
    return {"opm": np.ones(array.n, bool)}


@dataclass
class Sources:
    """Target sources (usable oct-6 vertices: inside the inner skull and >= 4 mm from its mesh) and
    the background grid (drawn from usable vertices)."""

    target: np.ndarray  # global full-resolution vertex indices
    grid: np.ndarray  # global indices of the 7-mm background grid
    grid_area: np.ndarray  # m^2 represented by each grid source (Voronoi cells)
    depth_mm: np.ndarray  # target depth below the dense scalp
    orientation_deg: np.ndarray  # angle to the local inner-skull normal (0 radial ... 90 tangential)
    lobe: np.ndarray
    region: np.ndarray  # Desikan-Killiany parcel of each target ('unknown' on the medial wall)
    dist_inner_skull_mm: np.ndarray


def make_sources(subject, cortex, rng) -> Sources:
    from . import plotting

    n_lh = int(np.sum(cortex.hemi == 0))
    oct_global = np.concatenate([subject.src[0]["vertno"], subject.src[1]["vertno"] + n_lh])
    target = oct_global[cortex.usable[oct_global]]
    valid = np.flatnonzero(cortex.usable)  # A-BEM-DIST: >= 4 mm from the inner-skull mesh
    grid = valid[goldenholz.poisson_disk(cortex.rr[valid], 0.007, rng)]
    grid_area = noisemodel.grid_areas(cortex.rr[grid], cortex.rr[valid], cortex.area[valid])
    depth = anatomy.depth_to_surface(cortex.rr[target], subject.scalp) * 1e3
    orient = anatomy.orientation_angle(cortex.rr[target], cortex.nn[target], subject.inner_skull)
    names = np.concatenate([plotting.read_freesurfer_annot(subject.labels / f"{h}.aparc.annot") for h in ("lh", "rh")])
    region = np.array([f"{'lh' if cortex.hemi[v] == 0 else 'rh'}.{names[v]}" for v in target], dtype=object)
    return Sources(target, grid, grid_area, depth, orient, plotting.lobe_of(names[target]), region,
                   cortex.dist_inner_skull[target] * 1e3)


def gains(array: Array, subject, cortex, points: np.ndarray, fullres_job: str | None = None,
          conductivity=BEM_CONDUCTIVITY, n_check: int = 12) -> np.ndarray:
    """Lead fields (n_channels, len(points)) for global vertex indices, from a full-resolution
    matrix if one exists for this array, else computed (cached) with the BEM. A full-resolution
    matrix is used only after ``n_check`` of its columns agree with a direct computation for
    this array (so a stale or mismatched matrix is never picked up silently)."""
    if fullres_job is not None and (paths.CACHE / "fullres" / f"{fullres_job}.npy").exists():
        full, col = fullres_matrix(array, subject, cortex, fullres_job, conductivity, n_check)
        if np.any(col[points] < 0):
            raise ValueError("points outside the valid full-resolution set")
        return np.asarray(full[:, col[points]], dtype=np.float64)
    return forward.chunked_discrete_gain(array.info, subject.trans, cortex.rr[points], cortex.nn[points],
                                         subject.bem_model(conductivity), coil_def=opm.coil_def_file(),
                                         label=array.name).astype(np.float64)


def fullres_matrix(array: Array, subject, cortex, job: str, conductivity=None, n_check: int = 12):
    """Memory-mapped full-resolution lead field (n_channels, n_valid) of ``array`` and the map from
    global vertex index to its column (-1 if not valid), after the fingerprint and column checks
    of ``opmsquid.fullres.load`` (the job fixes the BEM; ``conductivity`` is not used)."""
    from . import fullres  # imported here: fullres builds its arrays with this module

    return fullres.load(job, array.info, subject, cortex, n_check)


FULLRES_JOBS = {"squid": "neuromag_bem006", "opm_matched": "opm_bem006", "opm204": "opm204_bem006", "opm_dense": "opm_dense_bem006"}
ER_FILE = neuromag.ER_FILE


def measured_noise(squid_info: mne.Info, filt, bads, window=(-0.2, 0.0), trim_s: float = 2.0) -> dict:
    """Per-channel variance in the analysis band (SQUID channel order, native units^2, no SSP):
    'empty_room' (sample empty-room recording), 'baseline' (pre-stimulus windows of the sample
    task recording, filtered as one continuous record) and 'brain' = baseline - empty_room.
    Channels in ``bads`` are NaN. Also returns the empty-room environment fit."""
    names = squid_info.ch_names
    out = {}
    for key, fname in (("empty_room", ER_FILE), ("baseline", neuromag.RAW_FILE)):
        raw = mne.io.read_raw_fif(paths.SAMPLE_MEG / fname, preload=True, verbose=False)
        data = raw.get_data(picks=[raw.ch_names.index(n) for n in names])
        data = filt.apply(data - data.mean(axis=1, keepdims=True))
        sf, n_trim = raw.info["sfreq"], int(trim_s * raw.info["sfreq"])
        if key == "baseline":
            ev = mne.find_events(raw, stim_channel="STI 014", verbose=False)[:, 0] - raw.first_samp
            a, b = int(round(window[0] * sf)), int(round(window[1] * sf))
            segs = [data[:, e + a:e + b] for e in ev if e + a > n_trim and e + b < data.shape[1] - n_trim]
            seg = np.concatenate(segs, axis=1)
            out["n_windows"] = len(segs)
        else:
            seg = data[:, n_trim:-n_trim]
            er_raw = raw
        out[key] = np.mean(seg**2, axis=1)
        out[f"{key}_n_samples"] = seg.shape[1]
    bad = np.isin(names, list(bads))
    for key in ("empty_room", "baseline"):
        out[key][bad] = np.nan
    out["brain"] = out["baseline"] - out["empty_room"]
    out["environment"] = environment.fit_empty_room(er_raw, squid_info, filt, bads=list(bads))
    return out


def head_position_variants(info: mne.Info, subject, step: float = 0.005, pitch_deg: float = 5.0,
                           fit_clearance: float = 0.020, dewar_spacing: float = 0.018) -> dict:
    """Source-blind SQUID head positions (device-to-head transforms): the measured one; the head
    translated by +/-``step`` along each device axis; pitched by +/-``pitch_deg`` about the device
    x axis through the head origin; and 'well_fitted', moved up (device +z) until the scalp is
    ``fit_clearance`` from the nearest magnetometer coil. Each entry: (4x4 transform, minimum and
    median scalp-to-magnetometer distance [m], feasible = minimum >= ``dewar_spacing``)."""
    from scipy.spatial import cKDTree

    kinds = neuromag.channel_kinds(info)
    coil_dev = np.array([c["loc"][:3] for c, k in zip(info["chs"], kinds) if k == "mag"])
    mri_head = np.linalg.inv(subject.trans["trans"])
    tree = cKDTree(subject.scalp.rr @ mri_head[:3, :3].T + mri_head[:3, 3])
    base = info["dev_head_t"]["trans"]
    origin_dev = np.linalg.inv(base)[:3, 3]

    def moved(motion):  # rigid head motion (device frame): dev_head' = dev_head @ inv(motion)
        return base @ np.linalg.inv(motion)

    def translate(d):
        m = np.eye(4)
        m[:3, 3] = d
        return m

    def pitch(deg):
        c, s_ = np.cos(np.radians(deg)), np.sin(np.radians(deg))
        rot = np.eye(4)
        rot[1:3, 1:3] = [[c, -s_], [s_, c]]
        return translate(origin_dev) @ rot @ translate(-origin_dev)

    motions = {"measured": np.eye(4)}
    for ax, name in enumerate("xyz"):
        for sgn, lab in ((1, "+"), (-1, "-")):
            motions[f"{name}{lab}{step * 1e3:g}mm"] = translate(sgn * step * np.eye(3)[ax])
    motions[f"pitch+{pitch_deg:g}deg"] = pitch(pitch_deg)
    motions[f"pitch-{pitch_deg:g}deg"] = pitch(-pitch_deg)

    def dist(t):
        return tree.query(coil_dev @ t[:3, :3].T + t[:3, 3])[0]

    for dz in np.arange(0.0, 0.03, 0.0005):
        if dist(moved(translate([0, 0, dz]))).min() <= fit_clearance:
            motions["well_fitted"] = translate([0, 0, dz])
            break
    out = {}
    for name, m in motions.items():
        t = moved(m)
        d = dist(t)
        out[name] = dict(trans=t, min_dist=float(d.min()), median_dist=float(np.median(d)), feasible=bool(d.min() >= dewar_spacing),
                         motion=m)
    return out


def sample_covariance(cov: np.ndarray, n_samples: int, rng: np.random.Generator, shrink: bool = True) -> np.ndarray:
    """Covariance estimated from ``n_samples`` independent Gaussian draws with covariance ``cov``
    (Ledoit-Wolf shrinkage if ``shrink``, else the empirical estimate)."""
    ev, u = np.linalg.eigh(0.5 * (cov + cov.T))
    keep = ev > 1e-12 * ev.max()
    x = (u[:, keep] * np.sqrt(ev[keep])) @ rng.standard_normal((int(keep.sum()), n_samples))
    return metrics.ledoit_wolf_covariance(x)[0] if shrink else metrics.empirical_covariance(x)


def array_noise(array: Array, g_grid: np.ndarray, grid_area: np.ndarray, brain_scale: float, env: environment.EnvironmentModel,
                enbw: float, opm_asd: float, corr: tuple | None = None) -> noisemodel.ArrayNoise:
    """Noise model of one array: intrinsic white noise, calibrated cortical background (independent,
    or correlated with corr = (grid_points, length)), and the common room field."""
    if array.name == "squid":
        ivar = np.array([SQUID_ASD[k] ** 2 for k in array.kinds]) * enbw
    else:
        ivar = np.full(array.n, opm_asd**2 * enbw)
    mc = background.moment_covariance(grid_area) if corr is None else background.moment_covariance(grid_area, corr[0], corr[1])
    brain = brain_scale * background.sensor_covariance(g_grid, mc)
    basis = environment.external_basis(array.info, env.r0, array.coil_def)
    return noisemodel.ArrayNoise(ivar, brain, env.covariance(basis), basis)


def evaluate(topo: np.ndarray, noise: noisemodel.ArrayNoise, chans: np.ndarray, condition: str,
             cov_est: dict | None = None) -> dict:
    """Peak-channel SNR, mean-power SNR (dB) and known-topography detectability of topographies
    (n_channels, n_sources) restricted to channel set ``chans`` with the oracle covariance, plus
    the plug-in detectability for each estimated covariance in ``cov_est`` (label -> full-array
    covariance for this condition)."""
    s = noise.signal(topo, condition)[chans]
    c = noise.covariance(condition)[np.ix_(chans, chans)]
    var = np.diag(c)
    out = dict(peak=metrics.peak_channel_snr(s, var), meanpow_db=metrics.mean_power_snr_db(s, var),
               detect=metrics.detectability(s, c))
    for label, ce in (cov_est or {}).items():
        out[f"detect_{label}"] = metrics.plugin_detectability(s, c, ce[np.ix_(chans, chans)])
    return out
