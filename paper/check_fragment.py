#!/usr/bin/env python3
"""Render one LaTeX template fragment of the manuscript with the facts and compile it alone (a writer's check).

Renders paper/<fragment>.tex.j2 (template syntax of build_paper.py: << F.name >> for every number), wraps it in the
shared preamble with the supplementary numbering (S1, Fig. S1, Table S1), and compiles it with tectonic in
paper/build/check/. Reports unknown facts, LaTeX errors, glyphs missing from the font and undefined references other
than those to sections, figures and tables outside the fragment.

Usage: .venv/bin/python paper/check_fragment.py paper/supp/S3_noise.tex.j2 [--main]
  --main: number as the main text (1., Fig. 1) instead of the supplement.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_paper as bp  # noqa: E402
import jinja2  # noqa: E402

WRAP = r"""\documentclass[preprint,12pt,authoryear]{elsarticle}
\input{../../preamble.tex}
%(numbering)s
\begin{document}
\linenumbers
%(body)s
\bibliographystyle{elsarticle-harv}
\bibliography{../../references}
\end{document}
"""
SUPP_NUMBERING = r"""\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\theequation}{S\arabic{equation}}"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("fragment", type=Path)
    ap.add_argument("--main", action="store_true")
    args = ap.parse_args(argv)
    frag = args.fragment.resolve()
    facts = bp.report_facts.build_facts(bp.ROOT)
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(str(HERE)), undefined=jinja2.StrictUndefined,
                             variable_start_string="<<", variable_end_string=">>", block_start_string="<%",
                             block_end_string="%>", comment_start_string="<#", comment_end_string="#>",
                             keep_trailing_newline=True, autoescape=False)
    env.globals["trim"] = bp.trim
    env.filters["words"] = bp.words
    view = bp.FactView(facts)
    try:
        body = env.get_template(str(frag.relative_to(HERE))).render(F=view)
    except jinja2.UndefinedError as e:
        raise SystemExit(f"UNKNOWN FACT: {e}")
    out = HERE / "build" / "check"
    out.mkdir(parents=True, exist_ok=True)
    tex = out / f"{frag.name.replace('.tex.j2', '')}.tex"
    tex.write_text(WRAP % {"numbering": "" if args.main else SUPP_NUMBERING, "body": body})
    r = subprocess.run(["tectonic", "-X", "compile", "--keep-logs", tex.name], cwd=out, capture_output=True, text=True)
    log = (out / f"{tex.stem}.log").read_text(errors="replace") if (out / f"{tex.stem}.log").exists() else ""
    problems = []
    if r.returncode:
        problems.append("LaTeX failed:\n" + ((r.stderr or "") + (r.stdout or ""))[-2500:])
    problems += [f"missing glyph: {g}" for g in sorted(set(re.findall(r"Missing character: There is no (\S+)", log)))]
    undef = sorted(set(re.findall(r"Reference `([^']+)' on page", log)))
    cites = sorted(set(re.findall(r"Citation `([^']+)' on page", log)))
    problems += [f"undefined citation: {c}" for c in cites]
    print(f"{frag.name}: {len(view.used)} facts used; output {tex.with_suffix('.pdf').relative_to(bp.ROOT)}")
    if undef:
        print("references to labels outside this fragment (fine if they exist elsewhere): " + ", ".join(undef))
    if problems:
        print("PROBLEMS:\n" + "\n".join(problems))
        raise SystemExit(1)
    print("OK")


if __name__ == "__main__":
    main()
