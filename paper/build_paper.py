#!/usr/bin/env python3
"""Build the NeuroImage manuscript: fill every number from the facts of scripts/report_facts.py, then compile.

The sources are LaTeX templates in this folder (manuscript.tex.j2, supplementary.tex.j2) in which every number read from
an analysis output is written as << F.name >>, a fact of scripts/report_facts.py (value as printed, the unrounded value
and its source "<result file> :: <key path>"). Rendering fails on an unknown fact. Each rendered source is made
self-contained (the shared preamble inlined) and written twice: to build/, where it is compiled, and to this folder
(manuscript.tex, supplementary.tex: the submission sources, with figures read from figures/). The PDFs go to this
folder, the highlights to highlights.txt (the journal asks for them as a separate file), and build/numbers_used.tsv
lists every fact used with its value and source. The build fails when a highlight exceeds 85 characters, when there
are not 3 to 5 of them, when the abstract exceeds ABSTRACT_MAX_WORDS words, or when the body holds a bullet list. tectonic
keeps the .aux and .bbl in build/, which paper/build_html.py reads for the HTML version.

Template syntax (Jinja2 with LaTeX-safe delimiters): << F.name >> prints a fact; <% ... %> is a statement; <# ... #> a
comment; << trim("Figure_S5", 0) >> gives the includegraphics options showing one part of a stacked figure. Fact values are escaped for LaTeX and their Unicode exponents written as math.

Usage: .venv/bin/python paper/build_paper.py [--no-pdf] [--only manuscript|supplementary]
"""
from __future__ import annotations

import argparse
import json
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
ABSTRACT_MAX_WORDS = 400  # the Guide for Authors of NeuroImage asks for at most 250
SUPERSCRIPT = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")
LATEX_SPECIAL = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$", "&": r"\&", "#": r"\#", "_": r"\_",
                 "%": r"\%", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}


UNSIGNED_RANGE = re.compile(r"^(\d[\d,]*(?:\.\d+)?)(%?) to (\d[\d,]*(?:\.\d+)?)(%?)$")


def latex_value(text: str) -> str:
    """A fact's printed value as LaTeX, in the style sheet's conventions: an unsigned range 'a to b' becomes 'a–b' (a
    signed one keeps 'to'), a thin space precedes %, special characters are escaped, '× 10⁻⁵' is set as math and the
    rest is left as Unicode. The digits are never changed."""
    text = str(text)
    m = UNSIGNED_RANGE.match(text)
    if m:
        text = f"{m.group(1)}–{m.group(3)}" + ("%" if (m.group(2) or m.group(4)) else "")
    parts = re.split(r"(×\s*10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+|\b\d+[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]*[⁰⁴⁵⁶⁷⁸⁹⁻⁺][⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]*)", text)
    out = []
    for part in parts:
        m = re.fullmatch(r"×\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)", part)
        b = re.fullmatch(r"(\d+)([⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)", part)
        if m:
            out.append(r"$\times 10^{" + m.group(1).translate(SUPERSCRIPT) + "}$")
        elif b:  # a bare power such as 10⁻¹⁰ (the text font has no ⁰ or ⁻)
            out.append("$" + b.group(1) + "^{" + b.group(2).translate(SUPERSCRIPT) + "}$")
        else:
            out.append("".join(LATEX_SPECIAL.get(c, c) for c in part))
    return re.sub(r"(\d)\\%", r"\1\\,\\%", "".join(out))


class FactView:
    """Attribute access to the facts' LaTeX values; records every fact used."""

    def __init__(self, facts: dict):
        self._facts = facts
        self.used: dict[str, dict] = {}

    def dp(self, name: str, digits: int) -> str:
        """A numeric fact printed from its unrounded value with the given number of decimals (no double rounding)."""
        if name not in self._facts:
            raise jinja2.UndefinedError(f"unknown fact F.{name}")
        self.used[name] = self._facts[name]
        return latex_value(f"{float(self._facts[name]['raw']):.{digits}f}")

    def __getattr__(self, name: str) -> str:
        if name.startswith("_"):
            raise AttributeError(name)
        if name not in self._facts:
            raise jinja2.UndefinedError(f"unknown fact F.{name}")
        self.used[name] = self._facts[name]
        return latex_value(self._facts[name]["value"])


NUMBER_WORDS = {str(i): w for i, w in enumerate("zero one two three four five six seven eight nine".split())}


def words(value: str) -> str:
    """A count below ten in words, as the journal writes counts in prose ('8' -> 'eight'); other values unchanged."""
    return NUMBER_WORDS.get(str(value).strip(), value)


STACKS = HERE / "figure_stacks.json"  # where each part of a stacked figure sits (written by export_figures.py)


def trim(figure: str, part: int) -> str:
    """\\includegraphics options that show one part of a stacked figure alone, e.g. a figure too tall for one page
    shown as a float and its continuation: 'trim=left bottom right top,clip' in bp."""
    lay = json.loads(STACKS.read_text())[figure]
    p = lay["parts"][part]
    side = (lay["width"] - p["width"]) / 2
    bottom = lay["height"] - p["top"] - p["height"]
    return f"trim={side:.2f}bp {bottom:.2f}bp {side:.2f}bp {p['top']:.2f}bp,clip"


PREAMBLE_INPUT = r"\input{../preamble.tex}"
GRAPHICSPATH = r"\graphicspath{{../figures/}{../../figures/}}"


def flatten(text: str, figures: str) -> str:
    """The rendered source with the shared preamble inlined and figures read from `figures`."""
    preamble = (HERE / "preamble.tex").read_text().replace(GRAPHICSPATH, r"\graphicspath{{%s}}" % figures)
    if PREAMBLE_INPUT not in text:
        raise SystemExit("rendered source does not \\input the shared preamble")
    return text.replace(PREAMBLE_INPUT, preamble.strip())


def plain(tex: str) -> str:
    """Approximate plain text of a LaTeX fragment (for counting words and characters)."""
    tex = re.sub(r"\\(?:citep|citet|citealp|ref)\{[^}]*\}", "X", tex)
    tex = re.sub(r"\$\\times 10\^\{([^}]*)\}\$", r"× 10^\1", tex)
    tex = tex.replace("\\,", " ").replace("~", " ").replace("\\%", "%").replace("--", "–")
    tex = re.sub(r"\\[a-zA-Z]+\*?(\[[^]]*\])?", "", tex)
    tex = re.sub(r"[{}$]", "", tex)
    return re.sub(r"\s+", " ", tex).strip()


def between(text: str, tag: str) -> str:
    m = re.search(rf"%% {tag}-BEGIN\n(.*?)%% {tag}-END", text, re.S)
    if not m:
        raise SystemExit(f"markers %% {tag}-BEGIN/END not found")
    return m.group(1)


def check_manuscript(text: str) -> list[str]:
    """The journal's limits on the front matter, and no bullet list in the body."""
    problems = []
    items = [plain(i) for i in re.findall(r"\\item\s+(.*)", between(text, "HIGHLIGHTS"))]
    if not 3 <= len(items) <= 5:
        problems.append(f"{len(items)} highlights (3 to 5 required)")
    problems += [f"highlight of {len(i)} characters (85 at most): {i}" for i in items if len(i) > 85]
    abstract = plain(re.sub(r"\\(begin|end)\{abstract\}", "", between(text, "ABSTRACT")))
    n_words = len(abstract.split())
    if n_words > ABSTRACT_MAX_WORDS:
        problems.append(f"abstract of {n_words} words ({ABSTRACT_MAX_WORDS} at most)")
    body = text.split(r"\end{frontmatter}", 1)[-1]
    if re.search(r"\\begin\{(itemize|enumerate)\}", body):
        problems.append("bullet or numbered list in the body")
    (HERE / "highlights.txt").write_text("Highlights\n\n" + "\n".join(f"- {i}" for i in items) + "\n")
    print(f"highlights: {len(items)}, longest {max(map(len, items))} characters; abstract: {n_words} words")
    return problems


def stray_percent(text: str) -> list[str]:
    """A bare % in the body comments out the rest of its source line, which here is a whole paragraph."""
    body = text.split(r"\begin{document}", 1)[-1]
    return [f"unescaped % in the body: ...{line[max(0, m.start() - 60):m.start() + 1]}"
            for line in body.split("\n") for m in re.finditer(r"(?<![\\%])%", line) if line[:m.start()].strip()]


def check_abstract_in_pdf(text: str, pdf: Path) -> None:
    """The abstract's closing words must reach the compiled PDF (a stray comment or macro can swallow them)."""
    abstract = plain(re.sub(r"\\(begin|end)\{abstract\}", "", between(text, "ABSTRACT")))
    tail = " ".join(abstract.split()[-6:])
    r = subprocess.run(["pdftotext", "-f", "1", "-l", "3", str(pdf), "-"], capture_output=True, text=True)
    shown = re.sub(r"\s+", " ", r.stdout.replace("\u2019", "'"))
    if r.returncode or tail.replace("\u2019", "'") not in shown:
        raise SystemExit(f"{pdf.name}: the abstract's closing words are not in the PDF: '{tail}'")


def render(doc: str, facts: dict) -> tuple[Path, dict]:
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(str(HERE)), undefined=jinja2.StrictUndefined,
                             variable_start_string="<<", variable_end_string=">>", block_start_string="<%",
                             block_end_string="%>", comment_start_string="<#", comment_end_string="#>",
                             keep_trailing_newline=True, autoescape=False)
    env.globals["trim"] = trim
    env.filters["words"] = words
    view = FactView(facts)
    text = env.get_template(f"{doc}.tex.j2").render(F=view)
    problems = stray_percent(text)
    if doc == "manuscript":
        problems += check_manuscript(text)
    if problems:
        raise SystemExit(f"{doc}: " + "; ".join(problems))
    BUILD.mkdir(exist_ok=True)
    out = BUILD / f"{doc}.tex"
    out.write_text(flatten(text, "../figures/"))
    (HERE / f"{doc}.tex").write_text(flatten(text, "figures/"))
    return out, view.used


def compile_pdf(tex: Path) -> Path:
    shutil.copy2(HERE / "references.bib", BUILD / "references.bib")
    r = subprocess.run(["tectonic", "-X", "compile", "--keep-logs", "--keep-intermediates", tex.name], cwd=BUILD,
                       capture_output=True, text=True)
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
            if doc == "manuscript":
                check_abstract_in_pdf(tex.read_text(), pdf)
            print(f"compiled {pdf.relative_to(ROOT)}")
    lines = ["fact\tvalue\tsource"] + [f"{k}\t{v['value']}\t{v['source']}" for k, v in sorted(used_all.items())]
    (BUILD / "numbers_used.tsv").write_text("\n".join(lines) + "\n")
    print(f"{len(used_all)} facts used; build/numbers_used.tsv")


if __name__ == "__main__":
    main()
