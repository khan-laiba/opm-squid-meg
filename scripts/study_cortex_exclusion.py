#!/usr/bin/env python3
"""How much cortex the source rule leaves out in each head (A-BEM-DIST: no source within 4 mm of the inner skull).

Sources and the cortical background are drawn only from white-surface vertices inside the inner-skull BEM surface and at
least 4 mm from its mesh (opmsquid.anatomy.FullResCortex.usable), because the boundary-element lead fields are not
converged closer to it. In the adult this drops 8.7 % of the cortex. In heads whose skull lies close to the cortex the
rule removes more, and it removes the cortex nearest the scalp, which is the cortex nearest the OPMs. This script
reports, for every anatomy of the pediatric analyses (the adult, the two scaled adults, the three infant templates and
the three school-aged children), the share of the white-surface area that is not usable, split into vertices outside
the inner skull and vertices within 4 mm of it, over the whole surface and over the vertices within 8 and 10 mm of the
scalp used (the near-scalp thresholds of the children's MRI check). Nothing is simulated; the native surfaces are
loaded exactly as the pediatric analyses load them (scripts/study_children_qc.py's loader), and the scaled adults are
built as scripts/g3b_pediatric_helmet.py builds them (the adult scaled about its MRI origin, the 4-mm rule applied on
the scaled meshes; the 2-year size factor is the 24-month template's occipitofrontal circumference over the adult's).

Output: results/g3b/g3b_cortex_exclusion.json
Usage: .venv/bin/python scripts/study_cortex_exclusion.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import mne  # noqa: E402
import numpy as np  # noqa: E402

from opmsquid import anatomy, io, pediatric  # noqa: E402

OUT = ROOT / "results" / "g3b" / "g3b_cortex_exclusion.json"
KEYS = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
SCALED = ("school", "size2yr")
G3B = ROOT / "results" / "g3b" / "g3b_summary.json"
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
    loaded = {k: qc.load(k, cfg)[0] for k in KEYS if k not in SCALED}
    adult = loaded["adult"]
    factors = dict(school=float(cfg["anatomy"]["school_age_scale"]),
                   size2yr=pediatric.head_size(loaded["infant2yr"])["ofc_mm"] / pediatric.head_size(adult)["ofc_mm"])
    notes = json.loads(G3B.read_text())["anatomies"]
    for k, f in factors.items():  # the factors the pediatric analysis used (its notes print them to 4 decimals)
        if not notes[k]["scale_note"].endswith(f"{f:.4f}"):
            raise ValueError(f"{k}: scale factor {f:.6f} is not the one recorded in {G3B.name} ({notes[k]['scale_note']})")
        loaded[k] = anatomy.scaled(adult, f, f"sample_x{f:.4f}")
    out = dict(status="DERIVED (the cortex left out by the source rule A-BEM-DIST in each head of the pediatric analyses; "
                      "surfaces loaded or scaled as those analyses do; nothing simulated)",
               scale_factors=factors,
               rule=f"usable: inside the inner-skull BEM surface and at least {anatomy.MIN_BEM_DISTANCE * 1e3:g} mm from its "
                    "5120-triangle mesh (opmsquid.anatomy.FullResCortex.usable)",
               near_scalp="vertices whose distance to the scalp used (the dense MRI scalp of each head; point-to-surface, as the "
                          "children's MRI check measures it) is below the threshold",
               anatomies={})
    for key in KEYS:
        subject = loaded[key]
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
