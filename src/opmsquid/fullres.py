"""Full-resolution lead fields (every valid white-surface vertex) shared by G1B, G1C, G2 and G4.

Each job is an array and a BEM. A matrix is stored in ``cache/fullres/<job>.npy`` with a sidecar
``<job>.meta.json`` holding a fingerprint of what it was computed from: channel names, coil types,
coil positions and orientations, the head-to-MRI transform, the BEM conductivities and the valid
vertex set. ``load`` refuses a matrix whose fingerprint differs from the current array, and checks
a few columns against a direct forward computation; ``compute`` recomputes a job whose stored
fingerprint is missing or stale instead of silently reusing it.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import mne
import numpy as np

from . import forward, g2, io, neuromag, opm, paths  # io: the code commit is fixed at import

BEMS = {"hunold": (0.33, 0.0042, 0.33), "bem006": (0.3, 0.006, 0.3), "bem06": (0.3, 0.06, 0.3)}
JOBS = {
    "neuromag4pt_hunold": ("T3-4pt", "hunold"),  # Neuromag T3, 4-point coils (Hunold et al.)
    "opm_hunold": ("opm", "hunold"),  # matched OPM array, Hunold conductivities
    "neuromag_bem006": ("T3", "bem006"),  # Neuromag T3, MNE default 3-layer BEM
    "opm_bem006": ("opm", "bem006"),
    "neuromag_bem06": ("T3", "bem06"),  # Goldenholz et al. as printed
    "opm204_bem006": ("opm204", "bem006"),
    "opm_dense_bem006": ("opm_dense", "bem006"),
}


def directory() -> Path:
    return paths.CACHE / "fullres"


def array_info(kind: str, subject) -> mne.Info:
    dig = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
    if kind == "opm":
        return g2.matched_opm(subject, dig).info
    if kind in g2.DENSE:
        return g2.dense_opm(subject, dig, kind).info
    return neuromag.load_info(kind)


def fingerprint(info: mne.Info, subject, bem: list, valid_idx: np.ndarray, cortex) -> str:
    """Hash of every input of the matrix: the channels (names, coil types, positions, orientations,
    device-to-head), the coil definitions in use (MNE's standard file and the OPM coil), the
    head-to-MRI transform, the BEM (conductivities and surface geometry), the source vertices with
    their positions and normals, and the MNE version."""
    h = hashlib.sha256()
    h.update(mne.__version__.encode())
    h.update(forward.coil_definitions_digest(opm.coil_def_file()).encode())
    h.update("\n".join(info.ch_names).encode())
    for ch in info["chs"]:
        h.update(np.int64(ch["coil_type"]).tobytes())
        h.update(np.round(np.asarray(ch["loc"][:12], float), 9).tobytes())
    dev_head = info["dev_head_t"]["trans"] if info["dev_head_t"] is not None else np.eye(4)
    h.update(np.round(np.asarray(dev_head, float), 9).tobytes())
    h.update(np.round(np.asarray(subject.trans["trans"], float), 9).tobytes())
    for surf in bem:
        h.update(np.float64(surf["sigma"]).tobytes())
        h.update(np.round(np.asarray(surf["rr"], float), 9).tobytes())
        h.update(np.asarray(surf["tris"], np.int64).tobytes())
    h.update(np.asarray(valid_idx, np.int64).tobytes())
    h.update(np.round(np.asarray(cortex.rr[valid_idx], float), 9).tobytes())
    h.update(np.round(np.asarray(cortex.nn[valid_idx], float), 9).tobytes())
    return h.hexdigest()


def _meta_path(job):
    return directory() / f"{job}.meta.json"


def stored_fingerprint(job: str) -> str | None:
    p = _meta_path(job)
    return json.loads(p.read_text()).get("fingerprint") if p.exists() else None


def compute(job: str, subject, cortex, force: bool = False, log=print) -> Path:
    """Compute (or keep, if its fingerprint matches) the lead field of ``job``."""
    kind, bem_name = JOBS[job]
    out_dir = directory()
    out_dir.mkdir(parents=True, exist_ok=True)
    idx = np.flatnonzero(cortex.valid)
    np.save(out_dir / "valid_index.npy", idx)
    info = array_info(kind, subject)
    bem = subject.bem_model(BEMS[bem_name])
    fp = fingerprint(info, subject, bem, idx, cortex)
    target = out_dir / f"{job}.npy"
    if target.exists() and not force and stored_fingerprint(job) == fp:
        log(f"{job}: up to date, kept")
        return target
    t0 = time.time()
    gain = forward.chunked_discrete_gain(info, subject.trans, cortex.rr[idx], cortex.nn[idx], bem, coil_def=opm.coil_def_file(),
                                         label=job)
    np.save(target, gain)
    (out_dir / f"{job}.channels.txt").write_text("\n".join(info.ch_names))
    _meta_path(job).write_text(json.dumps(dict(job=job, fingerprint=fp, conductivity=BEMS[bem_name], shape=list(gain.shape),
                                               head_surface_triangles=int(len(bem[0]["tris"])), computed_at_commit=io.RUN_COMMIT),
                                          indent=1))
    log(f"{job}: {gain.shape} in {time.time() - t0:.0f} s -> {target}")
    return target


def load(job: str, info: mne.Info, subject, cortex, n_check: int = 6):
    """Memory-mapped lead field of ``job`` (n_channels, n_valid) and the map from global vertex
    index to column (-1 if not valid), after checking the fingerprint against ``info`` (the
    array the caller will use) and ``n_check`` columns against a direct computation. A matrix
    without a sidecar (older cache) is accepted only if its columns match, and then gets one."""
    kind, bem_name = JOBS[job]
    target = directory() / f"{job}.npy"
    full = np.load(target, mmap_mode="r")
    valid_idx = np.load(directory() / "valid_index.npy")
    bem = subject.bem_model(BEMS[bem_name])
    fp = fingerprint(info, subject, bem, valid_idx, cortex)
    stored = stored_fingerprint(job)
    if stored is not None and stored != fp:
        raise ValueError(f"{job}.npy was computed for a different array or BEM (fingerprint mismatch): recompute it")
    col = np.full(cortex.n, -1)
    col[valid_idx] = np.arange(len(valid_idx))
    chk = valid_idx[np.linspace(0, len(valid_idx) - 1, n_check).astype(int)]
    direct = forward.chunked_discrete_gain(info, subject.trans, cortex.rr[chk], cortex.nn[chk], bem, coil_def=opm.coil_def_file())
    have = np.asarray(full[:, col[chk]], dtype=np.float64)
    if have.shape != direct.shape or not np.allclose(have, direct, rtol=1e-4, atol=1e-4 * np.abs(direct).max()):
        raise ValueError(f"{job}.npy does not match the current array ({info.ch_names[0]}...): recompute it")
    if stored is None:
        _meta_path(job).write_text(json.dumps(dict(job=job, fingerprint=fp, conductivity=BEMS[bem_name], shape=list(full.shape),
                                                   head_surface_triangles=int(len(bem[0]["tris"])),
                                                   computed_at_commit="unknown (fingerprint added on a verified load)"), indent=1))
    return full, col
