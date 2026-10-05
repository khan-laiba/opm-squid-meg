#!/usr/bin/env python3
"""Build the NeuroImage manuscript: fill every number from the facts of scripts/report_facts.py, then compile.

The sources are LaTeX templates in this folder (manuscript.tex.j2, supplementary.tex.j2) in which every number read from
an analysis output is written as << F.name >>, a fact of scripts/report_facts.py (value as printed, the unrounded value
and its source "<result file> :: <key path>"). Rendering fails on an unknown fact. The rendered sources go to build/,
the PDFs to this folder, and build/numbers_used.tsv lists every fact used with its value and source.

Template syntax (Jinja2 with LaTeX-safe delimiters): << F.name >> prints a fact; <% ... %> is a statement; <# ... #> a
comment. Fact values are escaped for LaTeX and their Unicode exponents written as math.

Usage: .venv/bin/python paper/build_paper.py [--no-pdf] [--only manuscript|supplementary]
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import jinja2

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BUILD = HERE / "build"
sys.path.insert(0, str(ROOT / "scripts"))

import report_facts  # noqa: E402

DOCS = ("manuscript", "supplementary")
SUPERSCRIPT = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")
LATEX_SPECIAL = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$", "&": r"\&", "#": r"\#", "_": r"\_",
                 "%": r"\%", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}


def latex_value(text: str) -> str:
    """A fact's printed value as LaTeX: special characters escaped, '× 10⁻⁵' as math, the rest left as Unicode."""
    parts = re.split(r"(×\s*10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)", str(text))
    out = []
    for part in parts:
        m = re.fullmatch(r"×\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)", part)
        if m:
            out.append(r"$\times 10^{" + m.group(1).translate(SUPERSCRIPT).replace("-", "-") + "}$")
        else:
            out.append("".join(LATEX_SPECIAL.get(c, c) for c in part))
    return "".join(out)


class FactView:
    """Attribute access to the facts' LaTeX values; records every fact used."""

    def __init__(self, facts: dict):
        self._facts = facts
        self.used: dict[str, dict] = {}

    def __getattr__(self, name: str) -> str:
        if name.startswith("_"):
            raise AttributeError(name)
        if name not in self._facts:
            raise jinja2.UndefinedError(f"unknown fact F.{name}")
        self.used[name] = self._facts[name]
        return latex_value(self._facts[name]["value"])


def render(doc: str, facts: dict) -> tuple[Path, dict]:
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(str(HERE)), undefined=jinja2.StrictUndefined,
                             variable_start_string="<<", variable_end_string=">>", block_start_string="<%",
                             block_end_string="%>", comment_start_string="<#", comment_end_string="#>",
                             keep_trailing_newline=True, autoescape=False)
    view = FactView(facts)
    text = env.get_template(f"{doc}.tex.j2").render(F=view)
    BUILD.mkdir(exist_ok=True)
    out = BUILD / f"{doc}.tex"
    out.write_text(text)
    return out, view.used


def compile_pdf(tex: Path) -> Path:
    for name in ("references.bib", "elsarticle-num-names.bst"):
        src = HERE / name
        if src.exists():
            shutil.copy2(src, BUILD / name)
    r = subprocess.run(["tectonic", "-X", "compile", "--keep-logs", tex.name], cwd=BUILD, capture_output=True, text=True)
    log = (r.stdout or "") + (r.stderr or "")
    (BUILD / f"{tex.stem}.tectonic.txt").write_text(log)
    if r.returncode:
        raise SystemExit(f"tectonic failed for {tex.name}; see {BUILD / (tex.stem + '.tectonic.txt')}\n" + log[-3000:])
    missing = sorted(set(re.findall(r"Missing character: There is no (\S+)", (BUILD / f"{tex.stem}.log").read_text()
                                    if (BUILD / f"{tex.stem}.log").exists() else "")))
    if missing:
        raise SystemExit(f"{tex.name}: glyphs missing from the font: {missing}")
    pdf = HERE / f"{tex.stem}.pdf"
    shutil.copy2(BUILD / f"{tex.stem}.pdf", pdf)
    return pdf


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-pdf", action="store_true", help="render the LaTeX sources only")
    ap.add_argument("--only", choices=DOCS, help="build one document")
    args = ap.parse_args(argv)
    facts = report_facts.build_facts(ROOT)
    problems = report_facts.check(facts, ROOT)
    if problems:
        raise SystemExit("facts check failed:\n" + "\n".join(problems[:20]))
    used_all: dict[str, dict] = {}
    for doc in ([args.only] if args.only else DOCS):
        if not (HERE / f"{doc}.tex.j2").exists():
            print(f"{doc}.tex.j2 not found; skipped")
            continue
        tex, used = render(doc, facts)
        used_all.update(used)
        print(f"rendered {tex.relative_to(ROOT)} ({len(used)} facts)")
        if not args.no_pdf:
            pdf = compile_pdf(tex)
            print(f"compiled {pdf.relative_to(ROOT)}")
    lines = ["fact\tvalue\tsource"] + [f"{k}\t{v['value']}\t{v['source']}" for k, v in sorted(used_all.items())]
    (BUILD / "numbers_used.tsv").write_text("\n".join(lines) + "\n")
    print(f"{len(used_all)} facts used; build/numbers_used.tsv")


if __name__ == "__main__":
    main()
