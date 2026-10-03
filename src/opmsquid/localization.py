"""IED source localization on a bounded subset (G4, NEW).

Distributed estimates follow MNE-Python's minimum-norm conventions (make_inverse_operator /
apply_inverse), written out on whitened lead fields so that any array (OPM coils included) and any
discrete source grid can be used:

* whitening W from the (estimated) noise covariance, rank-aware (opmsquid.metrics.whitener);
* depth weighting of the source covariance, R_ii = (||W g_i||^2)^(-depth), depth = 0.8;
* R scaled so that trace(G_w R G_w^T) equals the whitened rank, lambda^2 = 1 / SNR^2 (SNR = 3);
* MNE kernel K = R G_w^T (G_w R G_w^T + lambda^2 I)^-1 (applied to whitened data);
* dSPM: each row divided by its noise SD, ||K_i|| (the whitened noise is white).
The equivalent current dipole uses MNE's ``fit_dipole`` directly (see scripts/g4_localization.py).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import metrics


@dataclass
class MNEInverse:
    w: metrics.Whitener
    kernel: np.ndarray  # (n_src, rank), applies to whitened data
    noise_sd: np.ndarray  # (n_src,) dSPM normalisation

    @staticmethod
    def make(gain: np.ndarray, noise_cov: np.ndarray, snr: float = 3.0, depth: float = 0.8) -> "MNEInverse":
        w = metrics.whitener(noise_cov)
        gw = w.apply(gain)
        r = np.sum(gw**2, axis=0) ** (-depth)
        r *= w.rank / np.sum(np.sum(gw**2, axis=0) * r)  # trace(G_w R G_w^T) = rank
        a = (gw * r[None, :]) @ gw.T + np.eye(w.rank) / snr**2
        k = (r[:, None] * gw.T) @ np.linalg.inv(a)
        return MNEInverse(w, k, np.linalg.norm(k, axis=1))

    def apply(self, data: np.ndarray, method: str = "dSPM") -> np.ndarray:
        est = self.kernel @ self.w.apply(data)
        return est / self.noise_sd[:, None] if method == "dSPM" else est


def channel_views(arrays: dict, subsets=()) -> dict:
    """Localization views: every array with all its channels and the Neuromag array ('squid')
    restricted to each sensor type in ``subsets`` ('mag', 'grad'; the views are named 'squid_mag',
    'squid_grad'): name -> (array name, channel indices). The physical arrays come first, in their
    order. A sensor type without channels is an error."""
    views = {name: (name, np.arange(a.n)) for name, a in arrays.items()}
    if "squid" in arrays:
        for k in subsets:
            idx = np.flatnonzero(np.asarray(arrays["squid"].kinds) == k)
            if not len(idx):
                raise ValueError(f"no Neuromag channels of type {k!r} (expected 'mag' or 'grad')")
            views[f"squid_{k}"] = ("squid", idx)
    return views


def peak_error(estimate: np.ndarray, src_rr: np.ndarray, true_rr: np.ndarray) -> tuple[float, int]:
    """Distance [m] from the true position to the source with the largest |estimate|."""
    i = int(np.argmax(np.abs(estimate)))
    return float(np.linalg.norm(src_rr[i] - true_rr)), i


def support_recovery(estimate: np.ndarray, inside: np.ndarray) -> float:
    """Share of the N strongest sources that belong to the simulated patch, N being the number of
    inverse-grid sources inside it (``inside``: boolean mask of the grid sources that are patch
    members, e.g. ``np.isin(grid, members)``); 1 = support recovered, chance N / n_src."""
    inside = np.asarray(inside, bool)
    n = int(inside.sum())
    if n == 0:
        return float("nan")
    top = np.argsort(np.abs(estimate))[::-1][:n]
    return float(np.mean(inside[top]))


def perturb_trans(trans: np.ndarray, rng: np.random.Generator, shift: float = 0.002, angle_deg: float = 2.0) -> np.ndarray:
    """A bounded coregistration error: random translation of length ``shift`` and rotation of
    ``angle_deg`` about a random axis through the head origin, applied to a 4x4 transform."""
    ax = rng.normal(size=3)
    ax /= np.linalg.norm(ax)
    th = np.radians(angle_deg)
    k = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    rot = np.eye(3) + np.sin(th) * k + (1 - np.cos(th)) * k @ k
    d = rng.normal(size=3)
    m = np.eye(4)
    m[:3, :3] = rot
    m[:3, 3] = shift * d / np.linalg.norm(d)
    return trans @ m
