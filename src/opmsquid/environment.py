"""Environmental interference through the actual sensor response (G2).

The residual room field near the head is described by its low-order external expansion: a
homogeneous field (3 components, T) and a linear, curl- and divergence-free gradient (5
components, T/m) about a reference point r0. ``external_basis`` evaluates these 8 fields with
MNE's own coil integration points, orientations and weights (the 'accurate' level used by the
forward model), so gradiometers see exactly their physical response (zero for a uniform field,
the baseline-weighted gradient otherwise) and finite OPM cells average the field over the cell.

Coefficient statistics are estimated once, on the Neuromag empty-room recording
(``fit_empty_room``), and then applied to every array through its own basis. Expressing the
fields in the head frame of the measured head position makes the fitted room field comparable
between the rigid SQUID helmet and a head-mounted OPM array (static head, no motion).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import mne
from mne.forward._make_forward import _prep_meg_channels  # MNE 1.13.2 internal, pinned

N_EXT = 8
GRADIENT_BASIS = np.array([  # symmetric, traceless 3x3 gradient tensors [dB_i/dx_j]
    [[1, 0, 0], [0, 0, 0], [0, 0, -1]],
    [[0, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[0, 1, 0], [1, 0, 0], [0, 0, 0]],
    [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
    [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
], dtype=float)


def coil_integration(info: mne.Info, coil_def=None) -> list[dict]:
    """Per-channel integration points (head frame), orientations and weights (MNE 'accurate')."""
    ctx = mne.use_coil_def(coil_def) if coil_def else _Null()
    with ctx:
        return _prep_meg_channels(info, accuracy="accurate", verbose=False)["defs"]


def external_basis(info: mne.Info, r0=(0.0, 0.0, 0.04), coil_def=None) -> np.ndarray:
    """(n_channels, 8) response of each channel to unit homogeneous-field components
    (Bx, By, Bz; 1 T) and unit gradient components (GRADIENT_BASIS; 1 T/m) about r0 [m]."""
    r0 = np.asarray(r0, float)
    rows = []
    for c in coil_integration(info, coil_def):
        w, cos, rel = c["w"], c["cosmag"], c["rmag"] - r0
        uniform = w @ cos  # (3,)
        grad = np.array([w @ np.sum((rel @ g.T) * cos, axis=1) for g in GRADIENT_BASIS])
        rows.append(np.concatenate([uniform, grad]))
    return np.array(rows)


@dataclass
class EnvironmentModel:
    coef_cov: np.ndarray  # (8, 8) covariance of the external coefficients in the analysis band
    coef_timecourses: np.ndarray  # (8, n_times) fitted coefficients (for time-domain simulation)
    sfreq: float
    explained_fraction: dict  # empty-room variance explained per channel type
    r0: np.ndarray

    def covariance(self, basis: np.ndarray) -> np.ndarray:
        return basis @ self.coef_cov @ basis.T


def fit_empty_room(raw: mne.io.BaseRaw, squid_info: mne.Info, analysis_filter, r0=(0.0, 0.0, 0.04),
                   type_sigma: dict | None = None, trim_s: float = 2.0, bads=()) -> EnvironmentModel:
    """Fit the 8 external coefficients to empty-room data (analysis band) by weighted least
    squares on all good MEG channels (weights 1/type_sigma, default: per-type median RMS), using
    the sensor geometry of ``squid_info`` (same channels, measured head position). Channels in
    ``bads`` are left out of the fit and of the explained fractions."""
    picks = mne.pick_types(raw.info, meg=True, exclude=[])
    names = [raw.ch_names[p] for p in picks]
    if names != [ch["ch_name"] for ch in squid_info["chs"]]:
        raise ValueError("empty-room channels must match the SQUID info channel order")
    good = ~np.isin(names, list(bads))
    data = raw.get_data(picks=picks)
    data = analysis_filter.apply(data - data.mean(axis=1, keepdims=True))
    n_trim = int(trim_s * raw.info["sfreq"])
    data = data[good, n_trim:-n_trim]
    basis = external_basis(squid_info, r0)[good]
    kinds = np.array(["grad" if ch["unit"] == mne.io.constants.FIFF.FIFF_UNIT_T_M else "mag"
                      for ch in squid_info["chs"]])[good]
    if type_sigma is None:
        rms = np.sqrt(np.mean(data**2, axis=1))
        type_sigma = {k: float(np.median(rms[kinds == k])) for k in ("mag", "grad")}
    wts = np.array([1.0 / type_sigma[k] for k in kinds])
    coef, *_ = np.linalg.lstsq(basis * wts[:, None], data * wts[:, None], rcond=None)
    resid = data - basis @ coef
    explained = {k: float(1.0 - np.sum(resid[kinds == k] ** 2) / np.sum(data[kinds == k] ** 2)) for k in ("mag", "grad")}
    return EnvironmentModel(np.cov(coef), coef, float(raw.info["sfreq"]), explained, np.asarray(r0, float))


class _Null:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False
