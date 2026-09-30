"""Common noise model for the realistic adult comparison (G2):

    y_m(t) = L_target,m q_target(t) + L_bg,m q_bg(t) + E_m e(t) + eps_m(t)

For each array m, in the analysis band of a common filter:
* intrinsic (eps): white noise per channel type (variance = ASD^2 x ENBW of the composite filter);
* brain background (L_bg q_bg): common cortical source process (opmsquid.background), projected
  through each array's own lead fields, scale s^2 calibrated once (on Neuromag gradiometers);
* environment (E e): common external-field coefficients (opmsquid.environment) through each
  array's own coil response.
Conditions: 'intrinsic', 'intrinsic+brain', 'intrinsic+brain+env', and 'projected' (external
subspace removed by a noise-weighted projection applied identically to signal and noise).
Nothing forces equal sensor RMS, equal SNR or a fixed noise ratio between arrays.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

CONDITIONS = ("intrinsic", "intrinsic+brain", "intrinsic+brain+env", "projected")


@dataclass
class ArrayNoise:
    intrinsic_var: np.ndarray  # (n,)
    brain_cov: np.ndarray  # (n, n)
    env_cov: np.ndarray  # (n, n)
    ext_basis: np.ndarray  # (n, 8)

    def covariance(self, condition: str) -> np.ndarray:
        c = np.diag(self.intrinsic_var).astype(float)
        if condition in ("intrinsic+brain", "intrinsic+brain+env", "projected"):
            c = c + self.brain_cov
        if condition in ("intrinsic+brain+env", "projected"):
            c = c + self.env_cov
        if condition == "projected":
            p = self.projector()
            c = p @ c @ p.T
        return c

    def projector(self) -> np.ndarray:
        """P = I - E (E^T W E)^-1 E^T W, W = diag(1/intrinsic_var): removes the external-field
        subspace (homogeneous field + linear gradients) with noise-weighted least squares,
        so channels in T and T/m are combined consistently."""
        e = self.ext_basis
        w = 1.0 / self.intrinsic_var
        m = e.T @ (w[:, None] * e)
        return np.eye(len(w)) - e @ np.linalg.solve(m, e.T * w[None, :])

    def signal(self, topographies: np.ndarray, condition: str) -> np.ndarray:
        return self.projector() @ topographies if condition == "projected" else topographies


def grid_areas(grid_points: np.ndarray, vertex_points: np.ndarray, vertex_areas: np.ndarray) -> np.ndarray:
    """Cortical area represented by each background-grid source: vertex areas summed over each
    grid point's Voronoi cell (nearest grid point, Euclidean)."""
    from scipy.spatial import cKDTree

    owner = cKDTree(grid_points).query(vertex_points)[1]
    return np.bincount(owner, weights=vertex_areas, minlength=len(grid_points))
