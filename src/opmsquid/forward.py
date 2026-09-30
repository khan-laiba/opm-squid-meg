"""Forward models with a content-addressed cache.

``fixed_gain`` returns the gain of sources fixed along the cortical normal (MNE cortical patch
statistics), ``discrete_gain`` the gain of dipoles at arbitrary points with given orientations.
Gains are float32 (relative precision ~1e-7). Results are cached in ``cache/fwd/`` under a SHA-1 of everything that determines them: channel
names, coil types, coil geometry (``loc``), ``dev_head_t``, the head-MRI transform, source
positions/orientations/vertex numbers, BEM surfaces and conductivities, extra coil definitions
and the MNE version. MNE 1.13.2 always uses its 'accurate' coil integration for MEG forwards.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import mne

from . import paths


def _hash_update(h, *items):
    for it in items:
        if isinstance(it, np.ndarray):
            h.update(np.ascontiguousarray(it, dtype=np.float64 if it.dtype.kind == "f" else it.dtype).tobytes())
        elif isinstance(it, (list, tuple)):
            _hash_update(h, *it)
        else:
            h.update(repr(it).encode())


def cache_key(info: mne.Info, trans, bem_surfaces: list, sources: tuple, coil_def: Path | None, kind: str) -> str:
    h = hashlib.sha1()
    _hash_update(h, kind, mne.__version__)
    meg = [ch for ch in info["chs"] if ch["kind"] == mne.io.constants.FIFF.FIFFV_MEG_CH]
    _hash_update(h, [ch["ch_name"] for ch in meg], [int(ch["coil_type"]) for ch in meg],
                 np.array([ch["loc"] for ch in meg]), np.asarray(info["dev_head_t"]["trans"]))
    _hash_update(h, np.asarray(trans["trans"]) if trans is not None else "no-trans")
    for s in bem_surfaces:
        _hash_update(h, int(s["id"]), float(s["sigma"]), s["rr"], s["tris"])
    _hash_update(h, *sources)
    _hash_update(h, Path(coil_def).read_text() if coil_def else "no-extra-coils")
    return h.hexdigest()[:16]


def _load(key: str):
    f = paths.CACHE / "fwd" / f"{key}.npz"
    if f.exists():
        with np.load(f) as z:
            return z["gain"], json.loads(str(z["meta"]))
    return None


def _save(key: str, gain: np.ndarray, meta: dict) -> None:
    f = paths.CACHE / "fwd" / f"{key}.npz"
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_suffix(".tmp.npz")
    np.savez(tmp, gain=gain, meta=json.dumps(meta))
    tmp.replace(f)


_SOLUTIONS: dict = {}


def _solution(bem_surfaces: list):
    """BEM solution, memoised in-process by surface ids, conductivities and geometry."""
    h = hashlib.sha1()
    for s in bem_surfaces:
        _hash_update(h, int(s["id"]), float(s["sigma"]), s["rr"], s["tris"])
    key = h.hexdigest()
    if key not in _SOLUTIONS:
        if len(_SOLUTIONS) >= 4:
            _SOLUTIONS.pop(next(iter(_SOLUTIONS)))
        _SOLUTIONS[key] = mne.make_bem_solution(bem_surfaces, verbose=False)
    return _SOLUTIONS[key]


def fixed_gain(info: mne.Info, trans, src: mne.SourceSpaces, bem_surfaces: list,
               coil_def: Path | None = None, use_cache: bool = True, use_cps: bool = True
               ) -> tuple[np.ndarray, dict]:
    """Gain [T or T/m per A m] of cortical sources fixed along the cortical normal
    (``surf_ori=True, force_fixed=True``; cortical patch statistics normals if ``use_cps``),
    shape (n_meg_channels, n_sources)."""
    sources = ([s["vertno"] for s in src], [s["rr"][s["vertno"]] for s in src], [s["nn"][s["vertno"]] for s in src])
    key = cache_key(info, trans, bem_surfaces, sources, coil_def, f"fixed-cps{int(use_cps)}")
    if use_cache and (hit := _load(key)) is not None:
        return hit
    ctx = mne.use_coil_def(coil_def) if coil_def else _null()
    with ctx:
        fwd = mne.make_forward_solution(info, trans, src, _solution(bem_surfaces), meg=True, eeg=False,
                                        mindist=0.0, verbose=False)
    fwd = mne.convert_forward_solution(fwd, surf_ori=True, force_fixed=True, use_cps=use_cps, verbose=False)
    if any(not np.array_equal(a["vertno"], b["vertno"]) for a, b in zip(fwd["src"], src)):
        raise RuntimeError("the forward model dropped sources (outside the inner skull?)")
    meta = dict(ch_names=fwd["info"]["ch_names"], n_sources=int(fwd["nsource"]), key=key)
    gain = np.asarray(fwd["sol"]["data"], dtype=np.float32)  # 1e-7 relative precision, half the storage
    if use_cache:
        _save(key, gain, meta)
    return gain, meta


def discrete_gain(info: mne.Info, trans, rr_mri: np.ndarray, nn_mri: np.ndarray, bem_surfaces: list,
                  coil_def: Path | None = None, use_cache: bool = True) -> tuple[np.ndarray, dict]:
    """Gain of dipoles at points ``rr_mri`` (MRI frame) with unit orientations ``nn_mri``,
    shape (n_meg_channels, n_points)."""
    rr_mri, nn_mri = np.asarray(rr_mri, float), np.asarray(nn_mri, float)
    key = cache_key(info, trans, bem_surfaces, (rr_mri, nn_mri), coil_def, "discrete")
    if use_cache and (hit := _load(key)) is not None:
        return hit
    src = mne.setup_volume_source_space(pos=dict(rr=rr_mri, nn=nn_mri), verbose=False)
    ctx = mne.use_coil_def(coil_def) if coil_def else _null()
    with ctx:
        fwd = mne.make_forward_solution(info, trans, src, _solution(bem_surfaces), meg=True, eeg=False,
                                        mindist=0.0, verbose=False)
    if fwd["nsource"] != len(rr_mri):
        raise RuntimeError("the forward model dropped sources (outside the inner skull?)")
    g = np.asarray(fwd["sol"]["data"]).reshape(len(fwd["info"]["ch_names"]), -1, 3)
    mri_head = np.linalg.inv(trans["trans"]) if trans is not None else np.eye(4)
    nn_head = nn_mri @ mri_head[:3, :3].T  # orientations in the head frame of the forward
    gain = np.einsum("cik,ik->ci", g, nn_head).astype(np.float32)
    meta = dict(ch_names=fwd["info"]["ch_names"], n_sources=len(rr_mri), key=key)
    if use_cache:
        _save(key, gain, meta)
    return gain, meta


def chunked_discrete_gain(info: mne.Info, trans, rr_mri: np.ndarray, nn_mri: np.ndarray, bem_surfaces: list,
                          coil_def: Path | None = None, chunk: int = 20000, label: str = "") -> np.ndarray:
    """``discrete_gain`` over many points in cached chunks (float32 result, (n_channels, n))."""
    import time

    out = None
    for start in range(0, len(rr_mri), chunk):
        t0 = time.time()
        sl = slice(start, min(start + chunk, len(rr_mri)))
        g, _ = discrete_gain(info, trans, rr_mri[sl], nn_mri[sl], bem_surfaces, coil_def)
        if out is None:
            out = np.empty((g.shape[0], len(rr_mri)), np.float32)
        out[:, sl] = g
        if label:
            print(f"  {label}: sources {sl.start}-{sl.stop} of {len(rr_mri)} ({time.time() - t0:.1f} s)", flush=True)
    return out


class _null:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False
