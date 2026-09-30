"""Hunold et al. (2016) adult spike simulations, MEG part (G1B).

Paper definitions kept (docs/literature/hunold2016.md, published pages):
* sources at every full-resolution white-surface vertex, fixed normal to the cortex (p. 1148);
* depth = Euclidean distance to the nearest *node* of the 5120-triangle scalp BEM mesh; orientation
  = angle between the source normal and the normal of the nearest node of the 5120-triangle
  inner-skull mesh, 0 deg radial ... 90 deg tangential (p. 1148);
* bins 20-60 mm x 0-90 deg in 5 mm x 10 deg (p. 1150);
* 600 nAm peak dipoles; ~20 mm^2 patches grown within +/-10 deg of the seed orientation until
  the area first exceeds 20 mm^2, uniform moment density (patch totals ~612-678 nAm, median
  622) (pp. 1149-1151);
* spike waveform peak 600 nAm, sharp component ~80 ms, complex ~200 ms (Fig. 2b, p. 1150);
* background: a random 10 % of the cortical nodes, independent Gaussian noise filtered into
  the five EEG bands, weighted 0.4-0.6, summed and normalised to +/-10 nAm; 6 s at 1 kHz;
  one fixed realization reused across sources (pp. 1150-1151, Fig. 6);
* no sensor noise (brain-noise-only reference, p. 1158);
* SNR = spike amplitude / background amplitude in the single channel (per sensor type) with the
  largest noise-free spike amplitude; background amplitude over the 1 s before onset =
  mean(+|hilbert|) - mean(-|hilbert|) = 2 mean |hilbert| (p. 1151); linear scale.

Unstated in the paper, chosen here and labelled (U-HU items in docs/provenance_register.md):
waveform shape (PCHIP through the digitised Fig. 2b points), EEG band edges and weights, peak
normalisation, onset time, Hilbert edge handling, the SNR numerator (three variants), and source
sampling (stratified to the paper's per-bin counts instead of the unreproducible "dipole traces").
"""
from __future__ import annotations

import heapq

import numpy as np
from scipy import signal
from scipy.interpolate import PchipInterpolator
from scipy.spatial import cKDTree

# digitised spike-wave complex (Appendix A of docs/literature/hunold2016.md): (ms, nAm)
SPIKE_POINTS = [(0, 0), (4, -30), (8, -66), (11, -67), (13, 0), (15, 130), (17, 310), (20, 466), (22, 560),
                (24, 600), (29, 525), (33, 338), (39, -10), (44, -132), (49, -215), (54, -253), (60, -230),
                (66, -147), (71, -66), (75, 24), (80, 111), (87, 149), (93, 138), (100, 121), (110, 80),
                (120, 20), (130, -40), (140, -79), (150, -91), (160, -81), (170, -60), (180, -28), (190, -8),
                (198, 0)]
DEPTH_EDGES_MM = np.arange(20.0, 60.1, 5.0)
ORIENT_EDGES_DEG = np.arange(0.0, 90.1, 10.0)
# sources per bin in the paper (Fig. 3c,d; Appendix B): rows 20-25 ... 55-60 mm, cols 0-10 ... 80-90 deg
PAPER_DIPOLE_COUNTS = np.array([
    [45, 33, 29, 30, 24, 19, 17, 16, 30], [47, 43, 52, 63, 48, 56, 82, 91, 130],
    [67, 49, 56, 51, 63, 64, 95, 111, 177], [56, 55, 59, 66, 62, 75, 72, 106, 143],
    [65, 66, 66, 63, 64, 78, 64, 57, 72], [32, 35, 40, 41, 47, 51, 51, 57, 62],
    [42, 37, 43, 43, 36, 23, 33, 21, 17], [25, 24, 23, 17, 21, 17, 21, 15, 32]])
# digitised bin-mean SNR classes (lower edge of 0.5-wide colour classes; Appendix C)
PAPER_BIN_CLASS = {
    ("dipole", "mag"): np.array([
        [1.0, 2.0, 2.5, 3.5, 4.0, 5.0, 4.5, 6.0, 5.5], [1.0, 2.0, 2.5, 2.5, 3.5, 4.0, 4.5, 4.5, 4.5],
        [1.5, 1.5, 2.0, 2.5, 3.0, 3.5, 3.5, 3.5, 3.5], [1.5, 1.5, 1.5, 2.0, 2.5, 2.5, 3.0, 2.5, 3.0],
        [1.0, 1.5, 1.5, 1.5, 2.0, 2.0, 2.0, 2.0, 2.0], [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2.0, 2.0, 2.0],
        [1.0, 1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.0, 2.0], [1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5, 1.5]]),
    ("dipole", "grad"): np.array([
        [1.5, 3.0, 4.5, 5.5, 7.5, 9.5, 8.5, 9.5, 9.5], [1.5, 3.0, 4.0, 5.0, 6.0, 7.0, 7.0, 7.5, 7.5],
        [1.5, 2.5, 3.0, 4.0, 5.0, 5.0, 6.0, 6.0, 6.0], [1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0, 4.0, 4.5],
        [1.5, 1.5, 2.0, 2.0, 2.5, 2.5, 3.0, 3.0, 3.5], [1.5, 1.5, 1.5, 2.0, 2.0, 2.0, 2.0, 2.5, 2.5],
        [1.5, 1.5, 1.5, 1.5, 2.0, 2.0, 2.0, 2.5, 2.5], [1.0, 1.5, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5, 1.5]]),
    ("patch", "mag"): np.array([  # rows 20-25 ... 50-55 mm; cols 10-20 ... 80-90 deg
        [1.0, 1.5, 3.0, 3.5, 3.5, 3.0, 4.5, 4.5], [1.5, 1.5, 2.5, 2.5, 3.5, 4.0, 4.0, 4.0],
        [1.5, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 3.5], [1.5, 1.5, 1.5, 2.0, 2.5, 2.5, 2.5, 2.5],
        [1.0, 1.5, 1.5, 2.0, 2.0, 2.0, 2.0, 2.0], [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5],
        [1.0, 1.0, 1.0, 1.5, 1.5, 2.0, 1.5, 2.0]]),
    ("patch", "grad"): np.array([
        [1.5, 2.0, 5.5, 5.5, 7.0, 5.0, 8.0, 8.0], [2.5, 2.5, 4.0, 4.5, 6.0, 6.0, 7.0, 7.0],
        [2.0, 2.0, 3.5, 4.0, 4.5, 5.0, 5.5, 5.5], [2.0, 2.0, 2.5, 3.0, 2.5, 3.0, 3.5, 4.0],
        [1.5, 1.5, 2.0, 2.5, 2.5, 2.5, 3.0, 3.0], [1.5, 1.5, 2.0, 2.0, 2.0, 2.0, 2.5, 2.5],
        [1.0, 1.5, 1.0, 2.0, 1.5, 2.5, 2.5, 2.5]]),
}
EEG_BANDS_HZ = ((0.5, 4.0), (4.0, 8.0), (8.0, 13.0), (13.0, 30.0), (30.0, 45.0))  # U-HU-bands
BAND_WEIGHTS = (0.6, 0.55, 0.5, 0.45, 0.4)  # within the paper's 0.4-0.6 range (U-HU-bands)


def spike_waveform(fs: float = 1000.0, peak: float = 600e-9) -> np.ndarray:
    """Spike-wave complex sampled at fs from onset (0) to 0.2 s, max = ``peak`` [A m]."""
    t_ms, v = np.array(SPIKE_POINTS, float).T
    t = np.arange(0.0, 0.2, 1.0 / fs) * 1e3
    w = PchipInterpolator(t_ms, v)(np.clip(t, 0, t_ms[-1]))
    w[t > t_ms[-1]] = 0.0
    return w / w.max() * peak


def bem_node_descriptors(points: np.ndarray, normals: np.ndarray, scalp_nodes: np.ndarray,
                         skull_nodes: np.ndarray, skull_normals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Depth [mm] (nearest scalp BEM node) and orientation [deg] (angle to the normal of the nearest
    inner-skull BEM node, folded to 0-90)."""
    depth = cKDTree(scalp_nodes).query(points)[0] * 1e3
    j = cKDTree(skull_nodes).query(points)[1]
    cos = np.abs(np.sum(normals * skull_normals[j], axis=1)) / np.linalg.norm(skull_normals[j], axis=1)
    return depth, np.degrees(np.arccos(np.clip(cos, 0.0, 1.0)))


def bin_index(depth, orient, depth_edges=DEPTH_EDGES_MM, orient_edges=ORIENT_EDGES_DEG):
    """(row, col) bin indices (left-closed bins; -1 outside the grid)."""
    r = np.searchsorted(depth_edges, depth, side="right") - 1
    c = np.searchsorted(orient_edges, orient, side="right") - 1
    c = np.where(np.isclose(orient, orient_edges[-1]), len(orient_edges) - 2, c)
    ok = (r >= 0) & (r < len(depth_edges) - 1) & (c >= 0) & (c < len(orient_edges) - 1)
    return np.where(ok, r, -1), np.where(ok, c, -1)


def sample_by_bins(depth, orient, candidates, counts, rng):
    """Stratified sample of ``counts[r, c]`` candidate vertices per bin (fewer if unavailable).
    Returns (indices, achieved counts)."""
    r, c = bin_index(depth[candidates], orient[candidates])
    picks, achieved = [], np.zeros_like(counts)
    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            pool = candidates[(r == i) & (c == j)]
            k = min(int(counts[i, j]), len(pool))
            achieved[i, j] = k
            if k:
                picks.append(rng.choice(pool, k, replace=False))
    return np.concatenate(picks), achieved


def grow_patch(seed: int, adjacency, orient: np.ndarray, area: np.ndarray, valid: np.ndarray,
               target_area: float = 20e-6, window: float = 10.0):
    """Grow a patch from ``seed`` in order of geodesic (mesh-edge) distance, admitting only valid
    vertices whose orientation lies within +/-window of the seed's (0-2*window if the seed is
    below window, 90-2*window...90 above 90-window), until the area first exceeds target_area.
    Returns the member indices, or None if the admissible region is smaller than the target."""
    o = orient[seed]
    lo, hi = (0.0, 2 * window) if o < window else ((90.0 - 2 * window, 90.0) if o > 90.0 - window else (o - window, o + window))
    members, total = [], 0.0
    dist = {seed: 0.0}
    heap = [(0.0, seed)]
    done = set()
    indptr, indices, data = adjacency.indptr, adjacency.indices, adjacency.data
    while heap:
        d, v = heapq.heappop(heap)
        if v in done:
            continue
        done.add(v)
        members.append(v)
        total += area[v]
        if total > target_area:
            return np.array(members)
        for k in range(indptr[v], indptr[v + 1]):
            u = indices[k]
            if u in done or not valid[u] or not (lo <= orient[u] <= hi):
                continue
            nd = d + data[k]
            if nd < dist.get(u, np.inf):
                dist[u] = nd
                heapq.heappush(heap, (nd, u))
    return None


def background_timecourses(n_sources: int, n_times: int, fs: float, rng: np.random.Generator,
                           peak: float = 10e-9, bands=EEG_BANDS_HZ, weights=BAND_WEIGHTS) -> np.ndarray:
    """Independent EEG-like background moments (n_sources, n_times) [A m]: Gaussian noise filtered
    into each band (Butterworth order 4, zero phase), scaled to unit RMS, weighted, summed and
    normalised so that max |s| = ``peak`` for every source."""
    x = np.zeros((n_sources, n_times))
    for (f1, f2), w in zip(bands, weights):
        sos = signal.butter(4, [f1, f2], btype="bandpass", fs=fs, output="sos")
        b = signal.sosfiltfilt(sos, rng.standard_normal((n_sources, n_times)), axis=1)
        x += w * b / b.std(axis=1, keepdims=True)
    return x / np.abs(x).max(axis=1, keepdims=True) * peak


def background_amplitude(trace: np.ndarray, window: slice, segment_only: bool = False) -> np.ndarray:
    """Hunold background amplitude per channel over ``window``: mean(+env) - mean(-env) with
    env = |hilbert|, i.e. 2 mean(env). The analytic signal is computed on the whole (demeaned)
    trace and cropped to the window (default) or on the window segment alone."""
    if segment_only:
        seg = trace[:, window]
        env = np.abs(signal.hilbert(seg - seg.mean(axis=1, keepdims=True), axis=1))
    else:
        env = np.abs(signal.hilbert(trace - trace.mean(axis=1, keepdims=True), axis=1))[:, window]
    return 2.0 * env.mean(axis=1)


def spike_snr(topographies: np.ndarray, waveform: np.ndarray, background: np.ndarray, onset: int,
              a_bg: np.ndarray) -> dict:
    """Hunold SNR for each topography (n_channels, n_sources) within one sensor type.

    Channel c* = the channel with the largest noise-free spike amplitude. Numerators:
    'peak' = noise-free max |spike|; 'p2p' = noise-free peak-to-peak; 'noisy_peak' = max |spike +
    background| in the 200-ms spike window. Denominator: a_bg[c*] (1 s before onset)."""
    wmax, wp2p = np.abs(waveform).max(), waveform.max() - waveform.min()
    c = np.argmax(np.abs(topographies), axis=0)
    g = topographies[c, np.arange(topographies.shape[1])]
    seg = background[c, onset:onset + len(waveform)]  # (n_sources, n_w)
    noisy = np.abs(g[:, None] * waveform[None, :] + seg).max(axis=1)
    return dict(channel=c, peak=np.abs(g) * wmax / a_bg[c], p2p=np.abs(g) * wp2p / a_bg[c], noisy_peak=noisy / a_bg[c])


def bin_means(values: np.ndarray, rows: np.ndarray, cols: np.ndarray, shape) -> tuple[np.ndarray, np.ndarray]:
    mean, count = np.full(shape, np.nan), np.zeros(shape, int)
    for i in range(shape[0]):
        for j in range(shape[1]):
            m = (rows == i) & (cols == j)
            count[i, j] = m.sum()
            if m.any():
                mean[i, j] = values[m].mean()
    return mean, count
