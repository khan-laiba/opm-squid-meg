"""Realistic adult OPM-Neuromag comparison (G2): arrays, sources, noise models, metrics.

Arrays (sample subject, measured head position in the Neuromag helmet):
  squid    Neuromag T3, 102 magnetometers + 204 planar gradiometers (comparators: 'mag', 'grad',
           'combined')
  opm99    matched-site OPM array (coverage control; 99 of 102 Neuromag sites feasible)
  opm204   204 sites spread evenly over the densest feasible array (channel-budget control vs the
           204 gradiometers)
  opm_dense  densest feasible single-axis OPM array under the 17-mm packing rule ("full system";
           216 sites on the sample head)
A 306-channel single-axis OPM array does not fit on this head (A-OPM-PACK); it is reported as
infeasible rather than simulated.
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
    arrays["opm99"] = matched_opm(subject, digitisation, scalp_gap)
    for name in DENSE:
        arrays[name] = dense_opm(subject, digitisation, name, scalp_gap)
    return arrays


# dense single-axis arrays: (farthest-point scalp spacing [m] of the densest feasible array,
# number of sites taken from it by farthest-point sampling (None = all), role)
DENSE = {"opm204": (0.015, 204, "channel-budget control vs 204 gradiometers"),
         "opm_dense": (0.015, None, "densest feasible single-axis array (full system)")}


def matched_opm(subject: anatomy.Subject, digitisation: mne.Info, scalp_gap: float = 0.0) -> Array:
    arr, rep = opm.matched_to_neuromag(neuromag.load_info("T3"), subject.trans, subject.scalp, digitisation,
                                       scalp_gap=scalp_gap)
    return _opm_array("opm99", arr, dict(role="matched-site coverage control", **_rep(rep)))


def dense_opm(subject: anatomy.Subject, digitisation: mne.Info, name: str, scalp_gap: float = 0.0) -> Array:
    """Densest feasible single-axis array (17-mm packing rule), or an evenly spread subset of it
    with a fixed channel count (farthest-point sampling of the sensing centres)."""
    spacing, n_sites, role = DENSE[name]
    skin = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    arr, rep = opm.dense_array(subject.scalp, subject.trans, digitisation, spacing, scalp_gap=scalp_gap, outer_skin=skin)
    if n_sites is not None:
        if n_sites > len(arr.pos):
            raise ValueError(f"{name}: only {len(arr.pos)} feasible sites")
        arr = opm.subset(arr, opm.farthest_point_subset(arr.pos, n_sites), f"dense OPM subset ({n_sites} sites)")
        rep = dict(rep, n_sites=n_sites, min_spacing_mm=float(opm.min_spacing(arr.pos).min() * 1e3),
                   median_spacing_mm=float(np.median(opm.min_spacing(arr.pos)) * 1e3), subset_of=rep["n_sites"])
    return _opm_array(name, arr, dict(role=role, **rep))


def _rep(rep):
    return dict(n_sites=rep["n_kept"], excluded_sites=rep["excluded_sites"], extra_shift_mm=rep["extra_shift_mm"],
                min_spacing_mm=float(np.min(rep["min_spacing_mm"])))


def _opm_array(name, arr, meta) -> Array:
    info = opm.make_info(arr)
    meta.update(channels=len(arr.pos), axes="1 (scalp normal)", standoff_mm=arr.standoff * 1e3,
                scalp_gap_mm=arr.scalp_gap * 1e3)
    return Array(name, info, np.array(["mag"] * len(arr.pos)), opm.coil_def_file(), meta)


def channel_sets(array: Array) -> dict[str, np.ndarray]:
    if array.name == "squid":
        return {"mag": array.kinds == "mag", "grad": array.kinds == "grad", "combined": np.ones(array.n, bool)}
    return {"opm": np.ones(array.n, bool)}


@dataclass
class Sources:
    """Target sources (oct-6 vertices inside the inner skull) and the background grid."""

    target: np.ndarray  # global full-resolution vertex indices
    grid: np.ndarray  # global indices of the 7-mm background grid
    grid_area: np.ndarray  # m^2 represented by each grid source (Voronoi cells)
    depth_mm: np.ndarray  # target depth below the dense scalp
    orientation_deg: np.ndarray  # angle to the local inner-skull normal (0 radial ... 90 tangential)
    lobe: np.ndarray


def make_sources(subject, cortex, rng) -> Sources:
    from . import plotting

    n_lh = int(np.sum(cortex.hemi == 0))
    oct_global = np.concatenate([subject.src[0]["vertno"], subject.src[1]["vertno"] + n_lh])
    target = oct_global[cortex.valid[oct_global]]
    valid = np.flatnonzero(cortex.valid)
    grid = valid[goldenholz.poisson_disk(cortex.rr[valid], 0.007, rng)]
    grid_area = noisemodel.grid_areas(cortex.rr[grid], cortex.rr[valid], cortex.area[valid])
    depth = anatomy.depth_to_surface(cortex.rr[target], subject.scalp) * 1e3
    orient = anatomy.orientation_angle(cortex.rr[target], cortex.nn[target], subject.inner_skull)
    names = np.concatenate([plotting.read_freesurfer_annot(paths.SUBJECTS_DIR / subject.name / "label" / f"{h}.aparc.annot")
                            for h in ("lh", "rh")])
    return Sources(target, grid, grid_area, depth, orient, plotting.lobe_of(names[target]))


def gains(array: Array, subject, cortex, points: np.ndarray, fullres_job: str | None = None) -> np.ndarray:
    """Lead fields (n_channels, len(points)) for global vertex indices, from a full-resolution
    matrix if available, else computed (cached) with the 3-layer BEM."""
    if fullres_job is not None:
        f = paths.CACHE / "fullres" / f"{fullres_job}.npy"
        if f.exists():
            valid_idx = np.load(paths.CACHE / "fullres" / "valid_index.npy")
            col = np.full(cortex.n, -1)
            col[valid_idx] = np.arange(len(valid_idx))
            return np.asarray(np.load(f, mmap_mode="r")[:, col[points]], dtype=np.float64)
    return forward.chunked_discrete_gain(array.info, subject.trans, cortex.rr[points], cortex.nn[points],
                                         subject.bem_model(BEM_CONDUCTIVITY), coil_def=opm.coil_def_file(),
                                         label=array.name).astype(np.float64)


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


def evaluate(topo: np.ndarray, noise: noisemodel.ArrayNoise, chans: np.ndarray, condition: str) -> dict:
    """Peak-channel SNR, mean-power SNR (dB) and known-topography detectability of topographies
    (n_channels, n_sources) restricted to channel set ``chans``, oracle covariance."""
    s = noise.signal(topo, condition)[chans]
    c = noise.covariance(condition)[np.ix_(chans, chans)]
    var = np.diag(c)
    return dict(peak=metrics.peak_channel_snr(s, var), meanpow_db=metrics.mean_power_snr_db(s, var),
                detect=metrics.detectability(s, c))
