"""Goldenholz et al. (2009) cortical SNR mapping, MEG part (G1C).

Eq. 1 (p. 1079): SNR = 10 log10[(a^2 / N) sum_k b_k^2 / s_k^2] -- a channel-averaged power SNR
with the 1/N factor, using per-channel noise variances only (no whitening).
Noise model (Eq. 3, p. 1080): s_k^2 = s_s^2 (A A^T)_kk for independent, cortex-normal noise
sources on a ~7 mm grid, with s_s^2 fitted to recorded data: for each sensor type the median over
channels of recorded_var_k / (A A^T)_kk, averaged with channel-count weights (p. 1080).
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree


def eq1_snr_db(topographies: np.ndarray, noise_var: np.ndarray) -> np.ndarray:
    """Goldenholz Eq. 1 for topographies (n_channels, n_sources) = a b_k (source amplitude
    times unit gain), with the 1/N factor over the given channels."""
    t = np.asarray(topographies, dtype=np.float64)
    ratio = np.einsum("ks,k->s", t**2, 1.0 / np.asarray(noise_var, float)) / t.shape[0]
    with np.errstate(divide="ignore"):
        return 10.0 * np.log10(ratio)


def calibrate_source_variance(recorded_var: np.ndarray, aat_diag: np.ndarray, kinds: np.ndarray) -> tuple[float, dict]:
    """Paper rule (p. 1080): per sensor type median_k(recorded_var_k / (A A^T)_kk), averaged with
    channel-count weights. Returns (s_s^2, per-type medians)."""
    per_type, num, den = {}, 0.0, 0
    for k in np.unique(kinds):
        m = kinds == k
        per_type[str(k)] = float(np.median(recorded_var[m] / aat_diag[m]))
        num += m.sum() * per_type[str(k)]
        den += m.sum()
    return num / den, per_type


def poisson_disk(points: np.ndarray, min_dist: float, rng: np.random.Generator) -> np.ndarray:
    """Greedy Poisson-disk subset: visit points in random order and keep a point if no kept point
    lies within ``min_dist`` (Euclidean)."""
    order = rng.permutation(len(points))
    kept = []
    grid = {}
    cell = min_dist
    for i in order:
        p = points[i]
        key = tuple(np.floor(p / cell).astype(int))
        ok = True
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                        if np.sum((points[j] - p) ** 2) < min_dist**2:
                            ok = False
                            break
                    if not ok:
                        break
                if not ok:
                    break
            if not ok:
                break
        if ok:
            kept.append(i)
            grid.setdefault(key, []).append(i)
    return np.sort(np.array(kept))


def geodesic_patches(adjacency, centroids: np.ndarray, radius: float, valid: np.ndarray, batch: int = 64) -> list:
    """Members of each patch: valid vertices within geodesic (mesh-edge Dijkstra) distance
    ``radius`` of each centroid."""
    out = []
    for start in range(0, len(centroids), batch):
        d = dijkstra(adjacency, directed=False, indices=centroids[start:start + batch], limit=radius)
        for row in d:
            m = np.flatnonzero(np.isfinite(row) & (row <= radius))
            out.append(m[valid[m]])
    return out


def nearest_indices(points: np.ndarray, targets: np.ndarray) -> np.ndarray:
    return cKDTree(points).query(targets)[1]
