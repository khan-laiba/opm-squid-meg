#!/usr/bin/env python3
"""Every number printed in the report (report/report.md), read from the committed result files.

build_facts() merges the per-milestone fact modules:
  scripts/report_facts_g12.py  G1A, G1B, G1C, G2 (adult) and the adult regions
  scripts/report_facts_g3.py   G3A, G3B (pediatric) and the regions by head
  scripts/report_facts_g4.py   G4 (spike detection, localization, motion)
  scripts/report_facts_methods.py  Methods details added in the revision (noise model, OPM axes, time-domain
                               model, bootstrap resamples, seeds, location-level test resolution, parameter provenance)
  scripts/report_facts_rev.py  the revision analyses (noise-model sensitivity, the noise model against the
                               measured Neuromag covariance, the children's MRI quality check)
  scripts/report_facts_writer.py  facts added by the writers of the revised report (the adult's primary ratios as D in dB,
                               the OPM averaging volume's effect from the register)
Each module exposes facts(root: Path) -> dict[name, fact]; a fact is
  {"value": text exactly as printed in the report, "raw": the unrounded number(s) or text,
   "source": "results/<file> :: <key path>" or "... :: derived: <how>"}.
Names are lower_snake_case, unique across modules; the report template uses them as F.<name>.

Usage: .venv/bin/python scripts/report_facts.py [--json OUT] [--check]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULES = ("report_facts_g12", "report_facts_g3", "report_facts_g4", "report_facts_methods", "report_facts_rev", "report_facts_confirm", "report_facts_writer")


def _load(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    if not path.exists():
        return None
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_facts(root: Path = ROOT, allow_missing: bool = False) -> dict:
    """All facts, name -> {"value", "raw", "source"}; raises on a duplicate name or a missing module."""
    out: dict = {}
    for name in MODULES:
        mod = _load(name)
        if mod is None:
            if allow_missing:
                continue
            raise FileNotFoundError(f"scripts/{name}.py is missing")
        for key, fact in mod.facts(root).items():
            if key in out:
                raise ValueError(f"fact {key!r} defined twice (second time in {name})")
            if not isinstance(fact, dict) or not {"value", "raw", "source"} <= fact.keys():
                raise ValueError(f"fact {key!r} in {name} lacks value/raw/source")
            out[key] = dict(fact, module=name)
    return out


def check(facts: dict, root: Path = ROOT) -> list[str]:
    """Problems: an empty value, or a source whose result file does not exist."""
    problems = []
    for key, f in facts.items():
        if str(f["value"]).strip() == "":
            problems.append(f"{key}: empty value")
        src = str(f["source"]).split("::")[0].strip()
        for part in src.replace(",", " ").split():
            if part.startswith(("results/", "configs/", "docs/")) and not (root / part).exists():
                problems.append(f"{key}: source file {part} not found")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--json", help="write every fact to this JSON file")
    ap.add_argument("--check", action="store_true", help="fail on an empty value or a missing source file")
    args = ap.parse_args(argv)
    facts = build_facts(ROOT)
    if args.json:
        Path(args.json).write_text(json.dumps(facts, indent=1, default=str))
    problems = check(facts) if args.check else []
    print(f"{len(facts)} facts from {len(MODULES)} modules" + (f"; {len(problems)} problems" if args.check else ""))
    if problems:
        print("\n".join(problems))
        sys.exit(1)


if __name__ == "__main__":
    main()
