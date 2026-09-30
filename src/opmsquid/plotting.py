"""Cortical maps on the inflated FreeSurfer surface (no nibabel needed).

Views are orthographic along +-x with z up; back faces are culled and triangles painted back to
front (painter's algorithm). Camera at -x shows anterior to the left (left lateral / right medial
views); camera at +x shows anterior to the right. These conventions were checked against
anatomical landmarks in the legacy realistic-head review.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection


def read_freesurfer_surface(fname) -> tuple[np.ndarray, np.ndarray]:
    """Vertices [mm] and triangles of a FreeSurfer binary triangle surface."""
    with open(fname, "rb") as fh:
        if fh.read(3) != b"\xff\xff\xfe":
            raise ValueError(f"not a FreeSurfer triangle surface: {fname}")
        fh.readline()
        fh.readline()
        n_vert, n_face = np.fromfile(fh, ">i4", 2)
        verts = np.fromfile(fh, ">f4", 3 * n_vert).reshape(-1, 3).astype(float)
        faces = np.fromfile(fh, ">i4", 3 * n_face).reshape(-1, 3)
    return verts, faces


def read_freesurfer_curv(fname) -> np.ndarray:
    """Per-vertex values of a FreeSurfer 'new format' curvature file (e.g. ?h.sulc; > 0 in sulci)."""
    with open(fname, "rb") as fh:
        if fh.read(3) != b"\xff\xff\xff":
            raise ValueError(f"not a FreeSurfer curvature file: {fname}")
        n_vert, _, _ = np.fromfile(fh, ">i4", 3)
        return np.fromfile(fh, ">f4", n_vert).astype(float)


def inflated_views(subjects_dir, subject, src):
    """Per hemisphere: (inflated vertices of the source-space vertices, triangles in source-space
    indices, gyral mask) for rendering source-space values."""
    surf_dir = Path(subjects_dir, subject, "surf")
    hemis = []
    for h, name in enumerate(("lh", "rh")):
        verts, _ = read_freesurfer_surface(surf_dir / f"{name}.inflated")
        sulc = read_freesurfer_curv(surf_dir / f"{name}.sulc")
        s = src[h]
        remap = np.full(s["np"], -1)
        remap[s["vertno"]] = np.arange(s["nuse"])
        tris = remap[s["use_tris"]]
        hemis.append((verts[s["vertno"]], tris[(tris >= 0).all(axis=1)], (sulc[s["vertno"]] < 0).astype(float)))
    return hemis


def render_view(ax, verts, tris, values, gyral, camera_x: float, cmap, norm, contour_level=None, shade=0.2):
    d = np.array([camera_x, 0.0, 0.0])
    screen = np.column_stack([verts[:, 1] * camera_x, verts[:, 2]])
    tri_n = np.cross(verts[tris[:, 1]] - verts[tris[:, 0]], verts[tris[:, 2]] - verts[tris[:, 0]])
    cos_view = tri_n @ d / np.linalg.norm(tri_n, axis=1)
    t = tris[cos_view > 0]
    t = t[np.argsort(verts[t].mean(axis=1) @ d)]
    vals = np.nanmean(values[t], axis=1)
    rgba = cmap(norm(vals))
    rgba[np.isnan(vals)] = (0.8, 0.8, 0.8, 1.0)
    rgba[:, :3] *= (1.0 - shade + shade * gyral[t].mean(axis=1))[:, None]
    ax.add_collection(PolyCollection(screen[t], facecolors=rgba, edgecolors="none", antialiased=False, rasterized=True))
    if contour_level is not None:
        ax.tricontour(mtri.Triangulation(screen[:, 0], screen[:, 1], tris[cos_view > 0.35]), np.nan_to_num(values),
                      levels=[contour_level], colors="k", linewidths=0.5)
    ax.set_xlim(screen[:, 0].min() - 2, screen[:, 0].max() + 2)
    ax.set_ylim(screen[:, 1].min() - 2, screen[:, 1].max() + 2)
    ax.set_aspect("equal")
    ax.axis("off")


def cortex_map_figure(hemis, rows: list[tuple[str, np.ndarray]], n_lh: int, cmap, norm, title: str, cbar_label: str,
                      fname, contour_level=None):
    """One row per (label, values over the source space), four views per row."""
    fig, axs = plt.subplots(len(rows), 4, figsize=(11, 1.9 * len(rows) + 0.8), squeeze=False)
    fig.subplots_adjust(left=0.16, right=0.9, top=1 - 0.6 / (1.9 * len(rows) + 0.8), bottom=0.03, wspace=0.02, hspace=0.05)
    for i, (label, vals) in enumerate(rows):
        for j, (h, cam) in enumerate(((0, -1.0), (0, 1.0), (1, -1.0), (1, 1.0))):
            verts, tris, gyral = hemis[h]
            v = vals[:n_lh] if h == 0 else vals[n_lh:]
            render_view(axs[i, j], verts, tris, v, gyral, cam, cmap, norm, contour_level)
        axs[i, 0].text(-0.05, 0.5, label, transform=axs[i, 0].transAxes, ha="right", va="center", fontsize=8)
    for j, t in enumerate(("left, lateral", "left, medial", "right, medial", "right, lateral")):
        axs[0, j].set_title(t, fontsize=8)
    cax = fig.add_axes([0.915, 0.15, 0.012, 0.7])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(cbar_label, fontsize=8)
    fig.suptitle(title, fontsize=9)
    fig.savefig(fname, dpi=160)
    plt.close(fig)
