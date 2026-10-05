#!/usr/bin/env python3
"""How much cortex the source rule leaves out in each head (A-BEM-DIST: no source within 4 mm of the inner skull).

Sources and the cortical background are drawn only from white-surface vertices inside the inner-skull BEM surface and at
least 4 mm from its mesh (opmsquid.anatomy.FullResCortex.usable), because the boundary-element lead fields are not
converged closer to it. In the adult this drops 8.7 % of the cortex. In heads whose skull lies close to the cortex the
rule removes more, and it removes the cortex nearest the scalp, which is the cortex nearest the OPMs. This script
reports, for the adult, the three infant templates and the three school-aged children (the native anatomies; the two
scaled adults are the adult scaled down), the share of the white-surface area that is not usable, split into vertices
outside the inner skull and vertices within 4 mm of it, over the whole surface and over the vertices within 8 and 10 mm
of the scalp used (the near-scalp thresholds of the children's MRI check). Nothing is simulated; the surfaces are
loaded exactly as the pediatric analyses load them (scripts/study_children_qc.py's loader).

Output: results/g3b/g3b_cortex_exclusion.json
Usage: .venv/bin/python scripts/study_cortex_exclusion.py
"""
from __future__ import annotations

import importlib.util
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import mne  # noqa: E402
import numpy as np  # noqa: E402

from opmsquid import anatomy, io  # noqa: E402

OUT = ROOT / "results" / "g3b" / "g3b_cortex_exclusion.json"
KEYS = ("adult", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
NEAR_MM = (8.0, 10.0)


def _qc():
    spec = importlib.util.spec_from_file_location("study_children_qc", ROOT / "scripts" / "study_children_qc.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["study_children_qc"] = mod
    spec.loader.exec_module(mod)
    return mod


def shares(area: np.ndarray, valid: np.ndarray, usable: np.ndarray, sel: np.ndarray) -> dict:
    a = area[sel]
    tot = float(a.sum())
    return dict(n_vertices=int(sel.sum()), area_cm2=tot * 1e4,
                usable_share=float(a[usable[sel]].sum() / tot) if tot else None,
                outside_inner_skull_share=float(a[~valid[sel]].sum() / tot) if tot else None,
                within_4mm_of_inner_skull_share=float(a[valid[sel] & ~usable[sel]].sum() / tot) if tot else None,
                n_usable_vertices=int(usable[sel].sum()))


def main():
    mne.set_log_level("WARNING")
    qc = _qc()
    cfg = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    out = dict(status="DERIVED (the cortex left out by the source rule A-BEM-DIST in each native head; surfaces loaded as the "
                      "pediatric analyses load them; nothing simulated)",
               rule=f"usable: inside the inner-skull BEM surface and at least {anatomy.MIN_BEM_DISTANCE * 1e3:g} mm from its "
                    "5120-triangle mesh (opmsquid.anatomy.FullResCortex.usable)",
               near_scalp="vertices whose distance to the scalp used (the dense MRI scalp of each head; point-to-surface, as the "
                          "children's MRI check measures it) is below the threshold",
               anatomies={})
    for key in KEYS:
        subject = qc.load(key, cfg)[0]
        cortex = anatomy.full_resolution(subject)
        valid, usable, area = cortex.valid, cortex.usable, cortex.area
        d_mm = qc.SignedDistance(subject.scalp.rr, subject.scalp.tris)(cortex.rr) * 1e3  # as the MRI check measures it
        rec = dict(subject=subject.name, whole_surface=shares(area, valid, usable, np.ones(cortex.n, bool)))
        for thr in NEAR_MM:
            rec[f"within_{thr:g}mm_of_scalp"] = shares(area, valid, usable, d_mm < thr)
        out["anatomies"][key] = rec
        w = rec["whole_surface"]
        print(f"{key:11s} usable {100 * w['usable_share']:.1f}% of {w['area_cm2']:.0f} cm2 (outside "
              f"{100 * w['outside_inner_skull_share']:.1f}%, within 4 mm {100 * w['within_4mm_of_inner_skull_share']:.1f}%); "
              f"within 8 mm of the scalp: {rec['within_8mm_of_scalp']['n_vertices']} vertices, "
              f"{rec['within_8mm_of_scalp']['n_usable_vertices']} usable", flush=True)
    io.write_json(out, OUT)
    print(f"-> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
