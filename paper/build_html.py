#!/usr/bin/env python3
"""Render the manuscript and its supplementary material as HTML pages (the HTML version of the paper on GitHub Pages).

The pages are converted from what the LaTeX build produced, without re-deciding anything that LaTeX and BibTeX
already decided:

  <build>/<doc>.tex   the rendered source, every number filled in (paper/build/manuscript.tex, supplementary.tex)
  <build>/<doc>.aux   the number of every \\label (sections, figures, tables, equations) and natbib's author-year
                      label of every cited key, as the PDF prints them
  <build>/<doc>.bbl   the reference list as BibTeX wrote it with elsarticle-harv (the PDF's reference list)
  <figures>/*.pdf     the figures, rasterized by pdftoppm at --dpi (default 200) with the \\includegraphics trim

and writes <out>/index.html (main text), <out>/supplementary.html, <out>/images/*.png and copies of the two PDFs
(the pages link to manuscript.pdf and supplementary.pdf beside them). Equations are typeset in the browser by
MathJax 3 from cdn.jsdelivr.net; nothing else is fetched, and the conversion itself needs no network. Every local
link and anchor of the output is checked (check_links); with --strict, a warning or a broken link is an error.

paper/build_paper.py keeps the .aux and .bbl in paper/build (tectonic --keep-intermediates). With --repo, the
script instead copies the rendered sources and references.bib of a repository into --work and compiles them there
with `tectonic -X compile --only-cached --keep-intermediates` (offline, from tectonic's cache); the repository is
only read.

Usage (from the repository root, after `.venv/bin/python paper/build_paper.py`):
  .venv/bin/python paper/build_html.py --build paper/build --figures paper/figures --out paper/build/html --strict
"""
from __future__ import annotations

import argparse
import html
import os
import re
import shutil
import struct
import subprocess
import time
import unicodedata
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from html.parser import HTMLParser
from pathlib import Path

MATHJAX_URL = "https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml.js"
TEXT_CM = 21.0 - 2 * 2.4  # \\textwidth of the A4 page with the preamble's margins
LANDSCAPE_CM = 29.7 - 2.4 - 2.6  # \\linewidth of a landscape page or sideways table
DOCS = {"manuscript": "index.html", "supplementary": "supplementary.html"}


# ---------------------------------------------------------------------------------------------------------- nodes
@dataclass
class Text:
    s: str


@dataclass
class Par:
    pass


@dataclass
class Cmd:
    name: str
    star: bool = False
    opts: list = field(default_factory=list)  # optional arguments: node lists, raw strings or None
    args: list = field(default_factory=list)  # mandatory arguments: node lists or raw strings


@dataclass
class Group:
    body: list


@dataclass
class Math:
    src: str
    display: bool = False
    env: str = ""


@dataclass
class Raw:
    html: str


@dataclass
class Env:
    name: str
    opts: list
    args: list
    body: list | None = None
    raw: str | None = None


# Argument signatures: s star, o optional [..] parsed, O optional raw, m mandatory parsed, v mandatory raw,
# p optional (..) raw. Commands not listed take no argument.
SIG = {
    "section": "som", "subsection": "som", "subsubsection": "som", "paragraph": "som",
    "label": "v", "ref": "v", "eqref": "v", "autoref": "v", "pageref": "v",
    "citep": "soov", "citet": "soov", "citealp": "soov", "citealt": "soov", "cite": "soov",
    "citeauthor": "soov", "citeyear": "soov", "citeyearpar": "soov",
    "textbf": "m", "textit": "m", "emph": "m", "texttt": "m", "textsc": "m", "textrm": "m", "textsf": "m",
    "textup": "m", "textnormal": "m", "textsuperscript": "m", "textsubscript": "m", "mbox": "m", "hbox": "m",
    "text": "m", "underline": "m", "url": "v", "href": "vm", "doi": "v", "path": "v", "nolinkurl": "v",
    "caption": "om", "includegraphics": "Ov", "captionsetup": "Ov", "multicolumn": "vvm", "cmidrule": "pv",
    "addlinespace": "O", "tablenotes": "m", "title": "om", "author": "Om", "address": "Om", "affiliation": "Om",
    "ead": "Om", "cortext": "Om", "corref": "v", "fntext": "Om", "fnref": "v", "journal": "m",
    "setlength": "vv", "renewcommand": "vv", "vspace": "sv", "hspace": "sv", "bibliographystyle": "v",
    "bibliography": "v", "item": "O", "\\": "sO", "bibinfo": "vm", "bibitem": "Ov", "natexlab": "m",
    "footnote": "Om", "thanks": "m", "rotatebox": "Ovm", "resizebox": "vvm", "scalebox": "vm", "phantom": "m",
}
ENV_SIG = {"figure": "O", "figure*": "O", "table": "O", "table*": "O", "sidewaystable": "O", "sidewaysfigure": "O",
           "tabular": "Ov", "tabular*": "vOv", "tabularx": "vv", "longtable": "Ov", "itemize": "O",
           "enumerate": "O", "minipage": "Ov", "thebibliography": "v"}
MATH_ENVS = {"equation", "equation*", "align", "align*", "gather", "gather*", "multline", "multline*",
             "eqnarray", "eqnarray*", "displaymath"}
RAW_ENVS = {"tabular", "tabular*", "tabularx", "longtable", "verbatim"}

CONTROL_SYMBOLS = {"%": "%", "&": "&", "#": "#", "_": "_", "$": "$", "{": "{", "}": "}", ",": "\u202f",
                   ";": "\u2005", ":": "\u2005", "!": "", " ": " ", "\n": " ", "\t": " ", "-": "", "@": "",
                   "/": ""}
WORD_SYMBOLS = {"ldots": "…", "dots": "…", "textendash": "–", "textemdash": "—", "textminus": "−", "S": "§",
                "dag": "†", "ddag": "‡", "textdagger": "†", "textdaggerdbl": "‡", "copyright": "©",
                "textbackslash": "\\", "textasciitilde": "~", "textasciicircum": "^", "quad": "\u2003",
                "qquad": "\u2003\u2003", "enspace": "\u2002", "thinspace": "\u202f", "nobreakspace": "\u00a0",
                "slash": "/", "hy": "-", "i": "ı", "o": "ø", "O": "Ø", "ss": "ß", "aa": "å", "AA": "Å",
                "ae": "æ", "AE": "Æ", "oe": "œ", "OE": "Œ", "l": "ł", "L": "Ł", "textdegree": "°", "textpm": "±",
                "texttimes": "×", "textbullet": "•", "textquoteright": "’", "textquoteleft": "‘",
                "textquotedblleft": "“", "textquotedblright": "”", "textmu": "µ", "checkmark": "✓",
                "newblock": " ", "DOIprefix": "", "URLprefix": "", "ArXivprefix": "", "sep": ", "}
ACCENTS = {"'": "\u0301", "`": "\u0300", "^": "\u0302", '"': "\u0308", "~": "\u0303", "=": "\u0304",
           ".": "\u0307", "u": "\u0306", "v": "\u030c", "H": "\u030b", "c": "\u0327", "k": "\u0328", "r": "\u030a"}
IGNORED = {"centering", "raggedright", "raggedleft", "flushleft", "footnotesize", "scriptsize", "tiny", "small",
           "normalsize", "large", "Large", "LARGE", "huge", "Huge", "FloatBarrier", "linenumbers",
           "nolinenumbers", "clearpage", "newpage", "pagebreak", "ContinuedFloat", "endfirsthead", "endhead",
           "endfoot", "endlastfoot", "toprule", "midrule", "bottomrule", "hline", "smallskip", "medskip",
           "bigskip", "protect", "relax", "arraybackslash", "noindent", "indent", "hfill", "vfill", "allowbreak",
           "nobreak", "null", "unskip", "ignorespaces", "LTcapwidth", "textwidth", "linewidth", "textheight",
           "baselineskip", "tabcolsep", "arraystretch", "thetable", "maketitle", "frenchspacing", "sloppy",
           "NoHyper", "endNoHyper", "appendix", "bibliographystyle", "label", "par", "strut", "nopagebreak"}
FONT_SWITCH = {"bfseries": "strong", "itshape": "em", "em": "em", "slshape": "em", "ttfamily": "code"}
SIZE_SWITCH = {"footnotesize", "scriptsize", "small", "tiny"}
SECTIONING = {"section": 1, "subsection": 2, "subsubsection": 3, "paragraph": 4}


class LatexError(ValueError):
    pass


def strip_comments(src: str) -> str:
    """Remove TeX comments: % to the end of the line, the newline and the next line's leading blanks."""
    return re.sub(r"(?<!\\)%[^\n]*(?:\n[ \t]*)?", "", src)


def read_group(s: str, i: int) -> tuple[str, int]:
    """The raw content of the balanced {...} starting at s[i] == '{', and the index after it."""
    if i >= len(s) or s[i] != "{":
        raise LatexError(f"expected {{ at {i}: {s[i:i + 40]!r}")
    depth, j = 0, i
    while j < len(s):
        c = s[j]
        if c == "\\":
            j += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    raise LatexError(f"unbalanced braces from {i}: {s[i:i + 60]!r}")


class Parser:
    """A small recursive-descent parser for the regular LaTeX that the build renders."""

    def __init__(self, src: str):
        self.s, self.i, self.n = src, 0, len(src)

    def fail(self, msg: str):
        raise LatexError(f"{msg} near {self.s[max(0, self.i - 60):self.i + 60]!r}")

    @classmethod
    def nodes(cls, src: str) -> list:
        return cls(src).parse()

    def parse(self, stop_env: str | None = None, stop_brace: bool = False) -> list:
        out: list = []
        buf: list[str] = []

        def flush():
            if buf:
                out.append(Text("".join(buf)))
                buf.clear()

        s = self.s
        while self.i < self.n:
            c = s[self.i]
            if c == "\\":
                if s.startswith("\\end{", self.i):
                    name, j = read_group(s, self.i + 4)
                    if name == stop_env:
                        flush()
                        self.i = j
                        return out
                    self.fail(f"unexpected \\end{{{name}}} inside {stop_env or 'top level'}")
                flush()
                node = self.command()
                if node is not None:
                    out.append(node)
            elif c == "{":
                flush()
                self.i += 1
                out.append(Group(self.parse(stop_brace=True)))
            elif c == "}":
                if not stop_brace:
                    self.fail("unbalanced }")
                flush()
                self.i += 1
                return out
            elif c == "$":
                flush()
                out.append(self.dollar_math())
            elif c == "~":
                buf.append("\u00a0")
                self.i += 1
            elif c in " \t\n":
                j, newlines = self.i, 0
                while j < self.n and s[j] in " \t\n":
                    newlines += s[j] == "\n"
                    j += 1
                self.i = j
                if newlines >= 2:
                    flush()
                    out.append(Par())
                else:
                    buf.append(" ")
            else:
                buf.append(c)
                self.i += 1
        if stop_env or stop_brace:
            self.fail(f"input ended inside {stop_env or 'a group'}")
        flush()
        return out

    # -- pieces
    def skip_blanks(self):
        """Spaces after a control word are not printed; one newline is a space, a blank line stays a paragraph."""
        s = self.s
        while self.i < self.n and s[self.i] in " \t":
            self.i += 1
        if self.i < self.n and s[self.i] == "\n":
            j = self.i + 1
            while j < self.n and s[j] in " \t":
                j += 1
            if j < self.n and s[j] == "\n":
                return
            self.i = j

    def optional(self) -> str | None:
        s = self.s
        if self.i >= self.n or s[self.i] != "[":
            return None
        depth, j = 0, self.i + 1
        while j < self.n:
            c = s[j]
            if c == "\\":
                j += 2
                continue
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
            elif c == "]" and depth == 0:
                raw = s[self.i + 1:j]
                self.i = j + 1
                return raw
            j += 1
        self.fail("unterminated [")

    def paren(self) -> str | None:
        s = self.s
        if self.i >= self.n or s[self.i] != "(":
            return None
        j = s.index(")", self.i)
        raw = s[self.i + 1:j]
        self.i = j + 1
        return raw

    def mandatory(self) -> str:
        s = self.s
        while self.i < self.n and s[self.i] in " \t\n":
            self.i += 1
        if self.i >= self.n:
            self.fail("missing argument")
        if s[self.i] == "{":
            raw, self.i = read_group(s, self.i)
            return raw
        if s[self.i] == "\\":
            m = re.compile(r"\\(?:[A-Za-z@]+|.)").match(s, self.i)
            self.i = m.end()
            return m.group(0)
        self.i += 1
        return s[self.i - 1]

    def until_end(self, name: str) -> str:
        end = "\\end{%s}" % name
        j = self.s.find(end, self.i)
        if j < 0:
            self.fail(f"no {end}")
        raw = self.s[self.i:j]
        self.i = j + len(end)
        return raw

    def dollar_math(self) -> Math:
        s = self.s
        if s.startswith("$$", self.i):
            j = s.find("$$", self.i + 2)
            raw, self.i = s[self.i + 2:j], j + 2
            return Math(raw, True)
        j = self.i + 1
        while j < self.n and not (s[j] == "$" and s[j - 1] != "\\"):
            j += 1
        if j >= self.n:
            self.fail("unterminated $")
        raw, self.i = s[self.i + 1:j], j + 1
        return Math(raw, False)

    def args(self, sig: str) -> tuple[bool, list, list]:
        star, opts, args = False, [], []
        for ch in sig:
            if ch == "s":
                if self.i < self.n and self.s[self.i] == "*":
                    star = True
                    self.i += 1
            elif ch in "oO":
                raw = self.optional()
                opts.append(raw if (raw is None or ch == "O") else Parser.nodes(raw))
            elif ch == "p":
                opts.append(self.paren())
            elif ch in "mv":
                raw = self.mandatory()
                args.append(raw if ch == "v" else Parser.nodes(raw))
        return star, opts, args

    def command(self):
        s = self.s
        j = self.i + 1
        if j >= self.n:
            self.fail("trailing backslash")
        if s[j].isalpha() or s[j] == "@":
            k = j
            while k < self.n and (s[k].isalpha() or s[k] == "@"):
                k += 1
            name, word = s[j:k], True
            self.i = k
        else:
            name, word = s[j], False
            self.i = j + 1
        if name == "begin":
            return self.environment()
        if name in ("(", "["):
            close = "\\)" if name == "(" else "\\]"
            k = s.find(close, self.i)
            raw, self.i = s[self.i:k], k + 2
            return Math(raw, name == "[")
        if name in ACCENTS and (not word or name in "uvHckr"):
            if word:
                self.skip_blanks()
            raw = self.mandatory()
            base = "".join(n.s if isinstance(n, Text) else WORD_SYMBOLS.get(getattr(n, "name", ""), "")
                           for n in Parser.nodes(raw))
            return Text(unicodedata.normalize("NFC", base[:1] + ACCENTS[name] + base[1:]))
        if not word and name in CONTROL_SYMBOLS:
            return Text(CONTROL_SYMBOLS[name])
        if word:
            self.skip_blanks()
        if name == "par":
            return Par()
        star, opts, args = self.args(SIG.get(name, ""))
        return Cmd(name, star, opts, args)

    def environment(self):
        name = self.mandatory()
        _, opts, args = self.args(ENV_SIG.get(name, ""))
        if name in MATH_ENVS:
            return Math(self.until_end(name), True, name)
        if name in RAW_ENVS:
            return Env(name, opts, args, raw=self.until_end(name))
        return Env(name, opts, args, body=self.parse(stop_env=name))


# ------------------------------------------------------------------------------------------------- aux and bbl
def latex_plain(raw: str) -> str:
    """Plain Unicode text of a short LaTeX fragment (a label number or a natbib author label)."""
    out = []
    for node in Parser.nodes(re.sub(r"\\nobreakspace\s*\{\}", "~", raw)):
        out.append(_plain_node(node))
    return re.sub(r"\s+", " ", "".join(out)).strip()


def _plain_node(node) -> str:
    if isinstance(node, Text):
        return ligatures(node.s)
    if isinstance(node, Group):
        return "".join(_plain_node(n) for n in node.body)
    if isinstance(node, Cmd):
        if node.name in WORD_SYMBOLS:
            return WORD_SYMBOLS[node.name]
        return "".join("".join(_plain_node(n) for n in a) if isinstance(a, list) else a for a in node.args)
    if isinstance(node, Math):
        return re.sub(r"\\[A-Za-z]+|[{}^_]", "", node.src)
    return ""


def read_aux(path: Path) -> tuple[dict, dict]:
    """\\newlabel numbers and \\bibcite labels (year, short author list) from a LaTeX .aux file."""
    labels, cites = {}, {}
    text = path.read_text(encoding="utf-8")
    for m in re.finditer(r"^\\newlabel\{([^}]*)\}", text, re.M):
        try:
            outer, _ = read_group(text, m.end())
            number, _ = read_group(outer, 0)
        except LatexError:
            continue
        labels[m.group(1)] = latex_plain(number)
    for m in re.finditer(r"^\\bibcite\{([^}]*)\}", text, re.M):
        outer, _ = read_group(text, m.end())
        parts, j = [], 0
        while j < len(outer) and len(parts) < 4:
            if outer[j] == "{":
                g, j = read_group(outer, j)
                parts.append(g)
            else:
                j += 1
        year, short = parts[1], parts[2]
        if short.startswith("{") and short.endswith("}"):
            short = short[1:-1]
        cites[m.group(1)] = {"year": year, "short": short}
    return labels, cites


def read_bbl(path: Path) -> list[tuple[str, list]]:
    """The entries of a .bbl file in order: (key, parsed content)."""
    text = strip_comments(path.read_text(encoding="utf-8"))
    start, end = text.find("\\bibitem"), text.rfind("\\end{thebibliography}")
    body = text[start:end]
    entries = []
    for chunk in re.split(r"(?=\\bibitem\b)", body):
        if not chunk.strip():
            continue
        p = Parser(chunk)
        p.i = len("\\bibitem")
        p.optional()
        key = p.mandatory()
        entries.append((key.strip(), Parser.nodes(chunk[p.i:])))
    return entries


def ligatures(s: str) -> str:
    """TeX's text ligatures: dashes and quotes."""
    return (s.replace("---", "\u2014").replace("--", "\u2013").replace("``", "\u201c").replace("''", "\u201d")
            .replace("`", "\u2018").replace("'", "\u2019"))


def esc(s: str) -> str:
    return html.escape(s, quote=False)


def attr(s: str) -> str:
    return html.escape(s, quote=True)


def html_id(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", label).strip("-")


def strip_tags(h: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", "", h))


def plain_sentence(h: str) -> str:
    """The first sentence of rendered HTML as plain text (for alt text): math simplified, links as numbers."""
    t = strip_tags(h.replace("\x00R:", "\x00N:"))
    t = re.sub(r"\\\((.*?)\\\)", lambda m: re.sub(r"\\(mathrm|symbf|mathbf|text)\b|[{}\\]", "", m.group(1)), t)
    return re.split(r"(?<=\.)\s(?=[A-Z(])", t, maxsplit=1)[0]


def autolink(h: str) -> str:
    """Bare URLs in text outside links become links (the reference list prints some URLs as plain text)."""
    out, inside_a = [], False
    for part in re.split(r"(<[^>]+>)", h):
        if part.startswith("<"):
            if re.match(r"<a\b", part):
                inside_a = True
            elif part.startswith("</a"):
                inside_a = False
            out.append(part)
        elif inside_a:
            out.append(part)
        else:
            out.append(re.sub(r"\b([Hh]ttps?://[^\s<]+?)(?=[.,;)]?(?:\s|$))",
                              lambda m: f"<a class=\"url\" href=\"{attr('h' + m.group(1)[1:])}\">{m.group(1)}</a>",
                              part))
    return "".join(out)


# ------------------------------------------------------------------------------------------------------ figures
def png_size(path: Path) -> tuple[int, int]:
    with open(path, "rb") as f:
        head = f.read(24)
    return struct.unpack(">II", head[16:24])


def pdf_page_size(pdf: Path) -> tuple[float, float]:
    out = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True, check=True).stdout
    m = re.search(r"Page size:\s+([\d.]+) x ([\d.]+) pts", out)
    return float(m.group(1)), float(m.group(2))


UNIT_BP = {"bp": 1.0, "pt": 72 / 72.27, "mm": 72 / 25.4, "cm": 72 / 2.54, "in": 72.0, "": 1.0}


def parse_trim(options: str | None) -> tuple[float, float, float, float] | None:
    """The trim=left bottom right top of \\includegraphics options, in bp (None without trim)."""
    if not options:
        return None
    m = re.search(r"trim\s*=\s*\{?([^,}]+)\}?", options)
    if not m:
        return None
    vals = []
    for tok in m.group(1).split():
        t = re.fullmatch(r"(-?[\d.]+)\s*([a-z]*)", tok)
        vals.append(float(t.group(1)) * UNIT_BP[t.group(2)])
    return tuple(vals)


@dataclass
class ImageJob:
    pdf: Path
    trim: tuple | None
    png: Path


def rasterize(job: ImageJob, dpi: int) -> None:
    """One figure (or the trimmed part of one) as a PNG at dpi, with pdftoppm."""
    cmd = ["pdftoppm", "-png", "-r", str(dpi), "-singlefile"]
    if job.trim:
        left, bottom, right, top = job.trim
        width, height = pdf_page_size(job.pdf)
        k = dpi / 72.0
        cmd += ["-x", str(round(left * k)), "-y", str(round(top * k)),
                "-W", str(round((width - left - right) * k)), "-H", str(round((height - top - bottom) * k))]
    stem = job.png.with_suffix("")
    subprocess.run(cmd + [str(job.pdf), str(stem)], check=True, capture_output=True)


# ------------------------------------------------------------------------------------------------------ renderer
class Renderer:
    def __init__(self, doc: str, tex: str, aux: tuple[dict, dict], bib: list, figures: Path, images: Path,
                 links: dict | None = None):
        self.doc = doc
        self.aux_labels, self.cites = aux
        self.bib = bib
        self.figures, self.images = figures, images
        self.links = links or {}  # (kind, number) -> URL in the other document
        self.warnings: list[str] = []
        self.labels: dict[str, tuple[str, str]] = {}  # label -> (number as computed, html id)
        self.counters = {"section": 0, "subsection": 0, "subsubsection": 0, "figure": 0, "table": 0, "equation": 0}
        self.appendix = False
        self.landscape = False
        self.n_cites = 0
        self.last_number = ""  # the number a stray \label refers to (LaTeX's last \refstepcounter)
        self.jobs: list[ImageJob] = []
        self.used_ids: set[str] = set()
        preamble, _, rest = tex.partition("\\begin{document}")
        self.body = strip_comments(rest.rpartition("\\end{document}")[0])
        self.prefix = {}
        for m in re.finditer(r"\\renewcommand\{\\the(section|figure|table|equation)\}\{([^\\]*)\\arabic", preamble):
            self.prefix[m.group(1)] = m.group(2)
        m = re.search(r"\\journal\{([^}]*)\}", preamble)
        self.journal = m.group(1) if m else ""
        links, self.links = self.links, {}
        self.ref_text = {key: re.sub(r"\s+", " ", strip_tags(self.inline(nodes))).strip() for key, nodes in bib}
        self.links = links
        self.title_text = ""

    # -- diagnostics
    def warn(self, msg: str):
        if msg not in self.warnings:
            self.warnings.append(msg)

    def new_id(self, base: str) -> str:
        hid, k = base, 2
        while hid in self.used_ids:
            hid, k = f"{base}-{k}", k + 1
        self.used_ids.add(hid)
        return hid

    def define(self, label: str, number: str) -> str:
        hid = self.new_id(html_id(label))
        self.labels[label] = (number, hid)
        return hid

    # -- numbers
    def number(self, kind: str) -> str:
        c = self.counters
        if kind == "section":
            if self.appendix:
                return "Appendix\u00a0" + chr(ord("A") + c["section"] - 1)
            return f"{self.prefix.get('section', '')}{c['section']}"
        if kind == "subsection":
            return f"{self.number('section')}.{c['subsection']}"
        if kind == "subsubsection":
            return f"{self.number('subsection')}.{c['subsubsection']}"
        return f"{self.prefix.get(kind, '')}{c[kind]}"

    # -- inline content
    def text(self, s: str) -> str:
        out = esc(ligatures(s))
        if self.links:
            out = self.link_other(out)
        return out

    TO_SUPP = re.compile(r"()\b(Sections?|Figs?\.|Tables?|Eqs?\.)(\u00a0| )(\(?)(S\d+(?:\.\d+)*)((?:[A-Z](?![a-z]))?)"
                         r"((?:(?:–|, | and )S\d+(?:\.\d+)*)?)")
    TO_MAIN = re.compile(r"(main[ -]text(?:\u2019s)?,? \(?)(Sections?|Figs?\.|Tables?|Eqs?\.)(\u00a0| )(\(?)"
                         r"(\d+(?:\.\d+)*)((?:[A-Z](?![a-z]))?)((?:(?:–|, | and )\d+(?:\.\d+)*)?)")
    KIND = {"Section": "sec", "Sections": "sec", "Fig.": "fig", "Figs.": "fig", "Table": "tab", "Tables": "tab",
            "Eq.": "eq", "Eqs.": "eq"}

    def link_other(self, h: str) -> str:
        """Links 'Section S4', 'Fig. S6A', 'Table S5', 'Sections S1–S10' in the main text, and 'main text, Section 2.5'
        or 'main text, Eq. (3)' in the supplement, to the other page when the number is one of its labels."""

        def a(kind, number):
            url = self.links.get((kind, number))
            return f"<a class=\"xref ext\" href=\"{attr(url)}\">{number}</a>" if url else number

        def repl(m):
            kind = self.KIND[m.group(2)]
            if (kind, m.group(5)) not in self.links:
                return m.group(0)
            out = f"{m.group(1)}{m.group(2)}{m.group(3)}{m.group(4)}{a(kind, m.group(5))}{m.group(6)}"
            tail = m.group(7)
            if tail:
                sep = re.match(r"–|, | and ", tail).group(0)
                out += sep + a(kind, tail[len(sep):])
            return out

        pattern = self.TO_SUPP if self.doc == "manuscript" else self.TO_MAIN
        return pattern.sub(repl, h)

    def math_inline(self, src: str) -> str:
        m = re.fullmatch(r"\s*\^\{?(?:\\mathrm\{)?([A-Za-z0-9,*+−-]+)\}?\}?\s*", src)
        if m:  # a superscript marker such as the a of a table note
            return f"<sup>{esc(m.group(1).replace('-', '−'))}</sup>"
        return f"<span class=\"math\">\\({esc(src.strip())}\\)</span>"

    def inline(self, nodes: list) -> str:
        parts = []
        for k, node in enumerate(nodes):
            if isinstance(node, Text):
                parts.append(self.text(node.s))
            elif isinstance(node, Par):
                parts.append(" ")
            elif isinstance(node, Raw):
                parts.append(node.html)
            elif isinstance(node, Group):
                parts.append(self.inline(node.body))
            elif isinstance(node, Math):
                if node.display:
                    parts.append(self.equation(node))
                else:
                    parts.append(self.math_inline(node.src))
            elif isinstance(node, Cmd):
                if node.name in FONT_SWITCH:
                    tag = FONT_SWITCH[node.name]
                    parts.append(f"<{tag}>{self.inline(nodes[k + 1:])}</{tag}>")
                    break
                parts.append(self.command(node))
            elif isinstance(node, Env):
                if node.body is not None:
                    parts.append(" ".join(self.blocks(node.body)))
                else:
                    self.warn(f"environment {node.name} in running text not rendered")
        return "".join(parts)

    def command(self, c: Cmd) -> str:
        n = c.name
        a = c.args
        wrap = {"textbf": "strong", "textit": "em", "emph": "em", "texttt": "code", "textsuperscript": "sup",
                "textsubscript": "sub", "underline": "u"}
        if n in wrap:
            return f"<{wrap[n]}>{self.inline(a[0])}</{wrap[n]}>"
        if n in ("mbox", "hbox", "text", "textrm", "textsf", "textup", "textnormal", "natexlab"):
            return self.inline(a[0])
        if n == "textsc":
            return f"<span class=\"sc\">{self.inline(a[0])}</span>"
        if n in ("url", "nolinkurl", "path"):
            u = re.sub(r"\\([%#_&$~])", r"\1", a[0])
            return f"<a class=\"url\" href=\"{attr(u)}\">{esc(u)}</a>"
        if n == "href":
            u = re.sub(r"\\([%#_&$~])", r"\1", a[0])
            return f"<a href=\"{attr(u)}\">{self.inline(a[1])}</a>"
        if n == "doi":  # the preamble's \doi: \url{https://doi.org/#1}, after an empty \DOIprefix
            u = "https://doi.org/" + a[0]
            return f"<a class=\"doi\" href=\"{attr(u)}\">{esc(u)}</a>"
        if n.startswith("cite"):
            return self.cite(c)
        if n in ("ref", "eqref", "autoref"):
            key = a[0].strip()
            return ("(" if n == "eqref" else "") + f"\x00R:{key}\x00" + (")" if n == "eqref" else "")
        if n in ("\\", "newline", "linebreak"):
            return "<br>"
        if n == "bibinfo":
            return f"<span class=\"bib-{attr(a[0])}\">{self.inline(a[1])}</span>"
        if n in WORD_SYMBOLS:
            return esc(WORD_SYMBOLS[n])
        if n in IGNORED or n in ("setlength", "renewcommand", "vspace", "hspace", "captionsetup", "addlinespace",
                                 "cmidrule", "phantom"):
            return ""
        if n == "footnote":
            self.warn("\\footnote rendered inline in brackets")
            return f" [{self.inline(a[0])}]"
        self.warn(f"unknown command \\{n} (arguments rendered as text)")
        return "".join(self.inline(x) for x in a if isinstance(x, list))

    def cite(self, c: Cmd) -> str:
        return f"<span class=\"citation\">{self._cite(c)}</span>"

    def _cite(self, c: Cmd) -> str:
        """natbib author-year citations as elsarticle sets them: (A et al., 2022, 2025; B, 2024)."""
        keys = [k.strip() for k in c.args[0].split(",") if k.strip()]
        self.n_cites += 1
        o = [x for x in c.opts if x is not None]
        pre = post = ""
        if len(o) == 1:
            post = self.inline(o[0]).strip()
        elif len(o) == 2:
            pre, post = self.inline(o[0]).strip(), self.inline(o[1]).strip()
        groups: list[tuple[str, list]] = []
        for key in keys:
            info = self.cites.get(key)
            if info is None:
                self.warn(f"citation {key} not in the .aux (shown as ?)")
                groups.append(("?", [(key, "?")]))
                continue
            short = self.inline(Parser.nodes(info["short"]))
            year = self.inline(Parser.nodes(info["year"]))
            if groups and groups[-1][0] == short:
                groups[-1][1].append((key, year))
            else:
                groups.append((short, [(key, year)]))

        def a(key, inner):
            tip = self.ref_text.get(key, "")
            return f"<a class=\"cite\" href=\"#ref-{attr(html_id(key))}\" title=\"{attr(tip)}\">{inner}</a>"

        def years(items):
            return ", ".join(a(k, y) for k, y in items)

        name = c.name
        if name == "citeauthor":
            return "; ".join(a(items[0][0], short) for short, items in groups)
        if name in ("citeyear", "citeyearpar"):
            inner = "; ".join(years(items) for _, items in groups)
            return f"({inner})" if name == "citeyearpar" else inner
        if name in ("citet", "citealt"):
            open_, close = (" (", ")") if name == "citet" else (" ", "")
            parts = [f"{a(items[0][0], short)}{open_}{years(items)}" + (f", {post}" if post else "") + close
                     for short, items in groups]
            return "; ".join(parts)
        parts = []
        for short, items in groups:
            first_key, first_year = items[0]
            first = a(first_key, f"{short}, {first_year}")
            rest = ", ".join(a(k, y) for k, y in items[1:])
            parts.append(first + (", " + rest if rest else ""))
        inner = "; ".join(parts)
        if pre:
            inner = f"{pre} {inner}"
        if post:
            inner = f"{inner}, {post}"
        return inner if name in ("citealp", "citealt") else f"({inner})"

    # -- blocks
    def blocks(self, nodes: list) -> list[str]:
        out: list[str] = []
        para: list = []

        def flush():
            h = self.inline(para).strip()
            if h:
                out.append(f"<p>{h}</p>")
            para.clear()

        k = 0
        while k < len(nodes):
            node = nodes[k]
            if isinstance(node, Par):
                flush()
            elif isinstance(node, Cmd) and node.name in SECTIONING:
                flush()
                label, k = self.next_label(nodes, k + 1)
                out.append(self.heading(node, label))
                continue
            elif isinstance(node, Math) and node.display:
                flush()
                out.append(self.equation(node))
            elif isinstance(node, Env):
                flush()
                out.extend(self.environment(node, out))
            elif isinstance(node, Cmd) and node.name == "tablenotes":
                flush()
                self.attach_notes(out, self.inline(node.args[0]))
            elif isinstance(node, Cmd) and node.name == "bibliography":
                flush()
                out.append(self.references())
            elif isinstance(node, Cmd) and node.name == "appendix":
                flush()
                self.appendix = True
                self.counters["section"] = self.counters["subsection"] = 0
            elif isinstance(node, Cmd) and node.name == "label":
                hid = self.define(node.args[0].strip(), self.last_number)
                para.append(Raw(f"<span id=\"{hid}\"></span>"))
            elif isinstance(node, Group) and self.is_block_group(node):
                flush()
                first = next((n for n in node.body if not (isinstance(n, Text) and not n.s.strip())), None)
                if isinstance(first, Cmd) and first.name in SIZE_SWITCH and out and "<!--NOTES-->" in out[-1]:
                    self.attach_notes(out, self.inline([n for n in node.body if not isinstance(n, Par)]).strip())
                else:
                    out.extend(self.blocks(node.body))
            elif isinstance(node, Cmd) and node.name in ("FloatBarrier", "linenumbers", "bibliographystyle",
                                                          "clearpage", "newpage", "vspace", "centering",
                                                          "setlength", "renewcommand", "captionsetup"):
                pass
            else:
                para.append(node)
            k += 1
        flush()
        return out

    @staticmethod
    def is_block_group(g: Group) -> bool:
        return any(isinstance(n, (Par, Env)) or (isinstance(n, Math) and n.display)
                   or (isinstance(n, Cmd) and n.name in SECTIONING) for n in g.body)

    @staticmethod
    def next_label(nodes: list, k: int) -> tuple[str | None, int]:
        j = k
        while j < len(nodes) and isinstance(nodes[j], Text) and not nodes[j].s.strip():
            j += 1
        if j < len(nodes) and isinstance(nodes[j], Cmd) and nodes[j].name == "label":
            return nodes[j].args[0].strip(), j + 1
        return None, k

    def attach_notes(self, out: list[str], notes: str):
        if out and "<!--NOTES-->" in out[-1]:
            out[-1] = out[-1].replace("<!--NOTES-->", f"<div class=\"tablenotes\">{notes}</div>")
        else:
            self.warn("table notes without a preceding table")
            out.append(f"<div class=\"tablenotes\">{notes}</div>")

    def heading(self, c: Cmd, label: str | None) -> str:
        level = SECTIONING[c.name]
        title = self.inline(c.args[0]).strip()
        if c.star or level > 3:
            number = ""
        else:
            kind = c.name
            self.counters[kind] += 1
            for lower in ("subsection", "subsubsection")[level - 1:]:
                self.counters[lower] = 0
            number = self.number(kind)
        self.last_number = number or self.last_number
        if label:
            hid = self.define(label, number)
        else:
            hid = self.new_id(("sec-" + number.replace("\u00a0", "-")) if number else
                              re.sub(r"[^a-z0-9]+", "-", strip_tags(title).lower()).strip("-")[:60])
        tag = f"h{level + 1}"
        num = f"<span class=\"secnum\">{esc(number)}.</span> " if number else ""
        return f"<{tag} id=\"{hid}\">{num}{title}</{tag}>"

    def equation(self, m: Math) -> str:
        src = m.src
        labels = re.findall(r"\\label\{([^}]*)\}", src)
        src = re.sub(r"\\label\{[^}]*\}", "", src).strip()
        numbered = not m.env.endswith("*") and m.env not in ("", "displaymath") and "\\nonumber" not in src
        if m.env in ("align", "align*", "gather", "gather*", "multline", "multline*", "eqnarray", "eqnarray*"):
            body = f"\\begin{{{m.env.rstrip('*')}*}}{src}\\end{{{m.env.rstrip('*')}*}}"
        else:
            body = f"\\[{src}\\]"
        if numbered:
            self.counters["equation"] += 1
            number = self.last_number = self.number("equation")
            hid = self.define(labels[0], number) if labels else self.new_id(f"eq-{number}")
            return (f"<div class=\"equation\" id=\"{hid}\"><div class=\"math display\">{esc(body)}</div>"
                    f"<div class=\"eqno\">({esc(number)})</div></div>")
        return f"<div class=\"equation\"><div class=\"math display\">{esc(body)}</div></div>"

    def environment(self, e: Env, out: list[str]) -> list[str]:
        name = e.name
        if name in ("figure", "figure*", "sidewaysfigure"):
            return [self.figure(e)]
        if name in ("table", "table*", "sidewaystable"):
            return [self.table(e, wide=name == "sidewaystable")]
        if name == "longtable":
            return [self.longtable(e)]
        if name in ("tabular", "tabular*", "tabularx"):
            head, body, _ = self.tabular(e)
            return [f"<div class=\"table-wrap\"><table>{head}{body}</table></div><!--NOTES-->"]
        if name == "frontmatter":
            return [self.frontmatter(e)]
        if name == "center":
            if self.doc == "supplementary" and not any(c > 0 for c in self.counters.values()):
                return [self.supplement_header(e)]
            return ["<div class=\"center\">" + "".join(self.blocks(e.body)) + "</div>"]
        if name in ("landscape", "NoHyper", "document", "minipage", "flushleft"):
            outer, self.landscape = self.landscape, self.landscape or name == "landscape"
            inner = self.blocks(e.body)
            self.landscape = outer
            return inner
        if name in ("itemize", "enumerate", "highlights"):
            tag = "ol" if name == "enumerate" else "ul"
            return [f"<{tag}>" + "".join(f"<li>{i}</li>" for i in self.items(e.body)) + f"</{tag}>"]
        if name == "abstract":
            return ["<section class=\"abstract\"><h2>Abstract</h2>" + "".join(self.blocks(e.body)) + "</section>"]
        self.warn(f"environment {name} rendered as plain blocks")
        return self.blocks(e.body or [])

    def items(self, body: list) -> list[str]:
        items, cur = [], None
        for node in body:
            if isinstance(node, Cmd) and node.name == "item":
                if cur is not None:
                    items.append(cur)
                cur = []
            elif cur is not None:
                cur.append(node)
        if cur is not None:
            items.append(cur)
        return [self.inline(i).strip() for i in items]

    # -- front matter
    def frontmatter(self, e: Env) -> str:
        title, authors, addresses, cortexts = "", [], {}, {}
        highlights, abstract, keywords = [], "", []
        for node in e.body:
            if isinstance(node, Cmd) and node.name == "title":
                title = self.inline(node.args[0]).strip()
            elif isinstance(node, Cmd) and node.name == "author":
                refs = [r.strip() for r in (node.opts[0] or "").split(",") if r.strip()]
                cor = [n.args[0] for n in node.args[0] if isinstance(n, Cmd) and n.name == "corref"]
                name = self.inline([n for n in node.args[0] if not (isinstance(n, Cmd) and n.name == "corref")])
                authors.append({"name": name.strip(), "refs": refs, "cor": cor, "email": None})
            elif isinstance(node, Cmd) and node.name == "ead":
                if authors:
                    authors[-1]["email"] = strip_tags(self.inline(node.args[0])).strip()
            elif isinstance(node, Cmd) and node.name == "cortext":
                cortexts[(node.opts[0] or "").strip()] = self.inline(node.args[0]).strip()
            elif isinstance(node, Cmd) and node.name in ("address", "affiliation"):
                addresses[(node.opts[0] or "").strip()] = self.inline(node.args[0]).strip()
            elif isinstance(node, Env) and node.name == "highlights":
                highlights = self.items(node.body)
            elif isinstance(node, Env) and node.name == "abstract":
                abstract = "".join(self.blocks(node.body))
            elif isinstance(node, Env) and node.name == "keyword":
                kw, cur = [], []
                for n in node.body:
                    if isinstance(n, Cmd) and n.name == "sep":
                        kw.append(cur)
                        cur = []
                    else:
                        cur.append(n)
                kw.append(cur)
                keywords = [self.inline(k).strip() for k in kw if self.inline(k).strip()]
        self.title_text = strip_tags(title)
        h = [f"<header class=\"front\"><h1 class=\"title\">{title}</h1>"]
        names = []
        for au in authors:
            marks = [esc(r) for r in au["refs"]]
            marks += [f"<a href=\"#{attr(html_id('fn-' + c))}\">∗</a>" for c in au["cor"]]
            sup = f"<sup>{','.join(marks)}</sup>" if marks else ""
            names.append(f"<span class=\"author\">{au['name']}{sup}</span>")
        h.append("<p class=\"authors\">" + ", ".join(names) + "</p>")
        if addresses:
            h.append("<ul class=\"affiliations\">" + "".join(
                f"<li id=\"aff-{attr(html_id(k))}\"><sup>{esc(k)}</sup>{v}</li>" for k, v in addresses.items())
                     + "</ul>")
        for key, text in cortexts.items():
            owners = [au for au in authors if key in au["cor"]]
            mail = ""
            for au in owners:
                if au["email"]:
                    mail += (f" Email address: <a href=\"mailto:{attr(au['email'])}\">{esc(au['email'])}</a> "
                             f"({strip_tags(au['name'])})")
            h.append(f"<p class=\"corresp\" id=\"{attr(html_id('fn-' + key))}\"><sup>∗</sup>{text}{mail}</p>")
        h.append("</header>")
        if highlights:
            h.append("<section class=\"highlights\" id=\"highlights\"><h2>Highlights</h2><ul>"
                     + "".join(f"<li>{i}</li>" for i in highlights) + "</ul></section>")
        if abstract:
            h.append(f"<section class=\"abstract\" id=\"abstract\"><h2>Abstract</h2>{abstract}")
            if keywords:
                h.append("<p class=\"keywords\"><span class=\"kw-label\">Keywords:</span> "
                         + ", ".join(f"<span class=\"kw\">{k}</span>" for k in keywords) + "</p>")
            h.append("</section>")
        return "".join(h)

    def supplement_header(self, e: Env) -> str:
        title, subtitle, rest = "", "", []
        for node in e.body:
            if isinstance(node, Group) and any(isinstance(n, Cmd) and n.name == "Large" for n in node.body):
                title = self.inline(node.body).strip()
            elif isinstance(node, Group) and any(isinstance(n, Cmd) and n.name == "large" for n in node.body):
                subtitle = self.inline(node.body).strip()
            elif not (isinstance(node, Cmd) and node.name == "\\"):
                rest.append(node)
        self.title_text = strip_tags(f"{title}: {subtitle}")
        return (f"<header class=\"front\"><p class=\"kicker\">{title}</p><h1 class=\"title\">{subtitle}</h1>"
                f"<p class=\"authors\">{self.inline(rest).strip()}</p></header>")

    # -- figures
    def figure(self, e: Env) -> str:
        images, caption, label, continued = [], None, None, False
        for node in e.body:
            if isinstance(node, Cmd) and node.name == "includegraphics":
                images.append(self.image(node.args[0].strip(), node.opts[0]))
            elif isinstance(node, Cmd) and node.name == "caption":
                caption = node.args[0]
            elif isinstance(node, Cmd) and node.name == "label":
                label = node.args[0].strip()
            elif isinstance(node, Cmd) and node.name == "ContinuedFloat":
                continued = True
        if not continued:
            self.counters["figure"] += 1
        number = self.last_number = self.number("figure")
        if label:
            hid = self.define(label, number)
        else:
            hid = self.new_id(f"fig-{number}" + ("-continued" if continued else ""))
        cap = self.inline(caption or []).strip()
        alt = f"Fig. {number}{' (continued)' if continued else ''}. {plain_sentence(cap)}"
        imgs = "".join(
            f"<a class=\"img-link\" href=\"{attr(src)}\"><img src=\"{attr(src)}\" width=\"{w}\" height=\"{h}\" "
            f"alt=\"{attr(alt)}\" loading=\"lazy\" decoding=\"async\"></a>" for src, w, h in images)
        return (f"<figure class=\"figure\" id=\"{hid}\">{imgs}<figcaption><span class=\"label\">"
                f"Fig.\u00a0{esc(number)}.</span> {cap}</figcaption></figure>")

    def image(self, name: str, options: str | None) -> tuple[str, int, int]:
        pdf = self.figures / (name if name.endswith(".pdf") else name + ".pdf")
        if not pdf.exists():
            raise SystemExit(f"figure not found: {pdf}")
        trim = parse_trim(options)
        same = [j for j in self.jobs if j.pdf == pdf]
        match = next((j for j in same if j.trim == trim), None)
        if match is None:
            stem = pdf.stem if trim is None else f"{pdf.stem}_part{len([j for j in same if j.trim]) + 1}"
            match = ImageJob(pdf, trim, self.images / f"{stem}.png")
            self.jobs.append(match)
        return f"images/{match.png.name}", 0, 0

    # -- tables
    def table(self, e: Env, wide: bool = False) -> str:
        caption, label, head, body, ncols = None, None, "", "", 0
        for node in e.body:
            if isinstance(node, Cmd) and node.name == "caption":
                caption = node.args[0]
            elif isinstance(node, Cmd) and node.name == "label":
                label = node.args[0].strip()
            elif isinstance(node, Env) and node.name in RAW_ENVS:
                head, body, ncols = self.tabular(node, LANDSCAPE_CM if (wide or self.landscape) else TEXT_CM)
        self.counters["table"] += 1
        number = self.last_number = self.number("table")
        hid = self.define(label, number) if label else self.new_id(f"tab-{number}")
        notes = [n for n in e.body if isinstance(n, Cmd) and n.name == "tablenotes"]
        notes_html = "".join(f"<div class=\"tablenotes\">{self.inline(n.args[0])}</div>" for n in notes)
        cls = "table wide" if (wide or self.landscape) else "table"
        cls += f" cols{min(ncols, 8)}"
        return (f"<figure class=\"{cls}\" id=\"{hid}\"><figcaption><span class=\"label\">Table\u00a0{esc(number)}."
                f"</span> {self.inline(caption or []).strip()}</figcaption><div class=\"table-wrap\"><table>{head}"
                f"{body}</table></div>{notes_html}<!--NOTES--></figure>")

    def longtable(self, e: Env) -> str:
        raw = e.raw
        caption, label = None, None
        m = re.match(r"\s*\\caption", raw)
        if m:
            p = Parser(raw)
            p.i = m.end()
            p.optional()
            caption = Parser.nodes(p.mandatory())
            raw = raw[p.i:]
            lm = re.match(r"\s*\\label\{([^}]*)\}", raw)
            if lm:
                label = lm.group(1)
                raw = raw[lm.end():]
            raw = re.sub(r"^\s*\\\\(\[[^\]]*\])?", "", raw)
        parts = re.split(r"\\(?:endfirsthead|endhead|endfoot|endlastfoot)\b", raw)
        first, body_raw = (parts[0], parts[-1]) if len(parts) > 1 else ("", raw)
        line_cm = LANDSCAPE_CM if self.landscape else TEXT_CM
        cols = parse_colspec(e.args[-1], line_cm)
        head_rows = table_rows(first)
        body_rows = table_rows(body_raw)
        self.counters["table"] += 1
        number = self.last_number = self.number("table")
        hid = self.define(label, number) if label else self.new_id(f"tab-{number}")
        head = self.rows_html(head_rows, cols, "th")
        body = self.rows_html(body_rows, cols, "td")
        cls = "table wide" if self.landscape else "table"
        cls += f" cols{min(len(cols), 8)}"
        return (f"<figure class=\"{cls}\" id=\"{hid}\"><figcaption><span class=\"label\">Table\u00a0{esc(number)}."
                f"</span> {self.inline(caption or []).strip()}</figcaption><div class=\"table-wrap\"><table>"
                f"{colgroup(cols, False, line_cm)}<thead>{head}</thead><tbody>{body}</tbody></table></div>"
                f"<!--NOTES--></figure>")

    def tabular(self, e: Env, line_cm: float = TEXT_CM) -> tuple[str, str, int]:
        cols = parse_colspec(e.args[-1], line_cm)
        rows = table_rows(e.raw)
        split = next((k for k, r in enumerate(rows) if k > 0 and "midrule" in r["rules"]), 0)
        head = self.rows_html(rows[:split], cols, "th")
        body = self.rows_html(rows[split:], cols, "td", skip_first_rule=bool(split))
        cg = colgroup(cols, e.name == "tabularx", line_cm)
        return cg + (f"<thead>{head}</thead>" if head else ""), f"<tbody>{body}</tbody>", len(cols)

    def rows_html(self, rows: list[dict], cols: list[tuple], cell: str, skip_first_rule: bool = False) -> str:
        out = []
        for k, row in enumerate(rows):
            classes = []
            if "midrule" in row["rules"] and not (skip_first_rule and k == 0):
                classes.append("midrule")
            cells, col = [], 0
            for span, align, raw in row["cells"]:
                align = align or (cols[col][0] if col < len(cols) else "l")
                cls = [{"l": "", "c": "c", "r": "r"}[align]]
                if any(a <= col + 1 and col + span <= b for a, b in row["cmid_below"]):
                    cls.append("cmid")
                c = " ".join(x for x in cls if x)
                attrs = (f" colspan=\"{span}\"" if span > 1 else "") + (f" class=\"{c}\"" if c else "")
                cells.append(f"<{cell}{attrs}>{self.inline(Parser.nodes(raw)).strip()}</{cell}>")
                col += span
            if len(row["cells"]) == 1 and row["cells"][0][0] == len(cols) and len(cols) > 1:
                classes.append("group")
            cls = f" class=\"{' '.join(classes)}\"" if classes else ""
            out.append(f"<tr{cls}>{''.join(cells)}</tr>")
        return "".join(out)

    # -- references
    def references(self) -> str:
        links, self.links = self.links, {}
        items = []
        for key, nodes in self.bib:
            entry = autolink(self.inline(nodes).strip())
            items.append(f"<li id=\"ref-{attr(html_id(key))}\">{entry}</li>")
        self.links = links
        return ("<section class=\"references\" id=\"references\"><h2>References</h2><ul class=\"reflist\">"
                + "".join(items) + "</ul></section>")

    # -- the page
    def render(self) -> str:
        nodes = Parser.nodes(self.body)
        blocks = self.blocks(nodes)
        page = "\n".join(b.replace("<!--NOTES-->", "") for b in blocks)

        def ref(m):
            key = m.group(1)
            number = self.aux_labels.get(key)
            mine = self.labels.get(key)
            if mine is None:
                self.warn(f"\\ref{{{key}}}: label not found in the document")
                return "??"
            if number is None:
                self.warn(f"\\ref{{{key}}}: not in the .aux; using the computed number {mine[0]}")
                number = mine[0]
            elif number != mine[0].replace("\u00a0", " "):
                self.warn(f"\\ref{{{key}}}: computed number {mine[0]!r} differs from the .aux's {number!r}")
            return f"<a class=\"xref\" href=\"#{attr(mine[1])}\">{esc(number.replace(' ', chr(160)))}</a>"

        page = re.sub(r"\x00R:([^\x00]*)\x00", ref, page)
        page = re.sub(r"\x00N:([^\x00]*)\x00", lambda m: self.aux_labels.get(m.group(1)) or
                      self.labels.get(m.group(1), ("??",))[0], page)
        for key, number in self.aux_labels.items():
            mine = self.labels.get(key)
            if mine and mine[0] and number != mine[0].replace("\u00a0", " "):
                self.warn(f"label {key}: computed number {mine[0]!r} differs from the .aux's {number!r}")
        return page


CM = {"cm": 1.0, "mm": 0.1, "in": 2.54, "pt": 2.54 / 72.27, "bp": 2.54 / 72}


def width_cm(w: str, line_cm: float) -> float | None:
    """A column width such as 3.0cm or 0.44\\textwidth in cm (line_cm: the width of \\textwidth/\\linewidth)."""
    m = re.fullmatch(r"\s*([\d.]+)\s*(cm|mm|in|pt|bp)\s*", w)
    if m:
        return float(m.group(1)) * CM[m.group(2)]
    m = re.fullmatch(r"\s*([\d.]*)\s*\\(textwidth|linewidth|columnwidth|hsize)\s*", w)
    if m:
        return float(m.group(1) or 1) * line_cm
    return None


def parse_colspec(spec: str, line_cm: float = TEXT_CM) -> list[tuple[str, float | None]]:
    """(alignment 'l'/'c'/'r', width in cm or None) of each column of a tabular column specification, with this
    preamble's L{w} (ragged right), R{w} (ragged left) and Y (ragged-right X) column types."""
    cols, pending, k = [], None, 0
    while k < len(spec):
        c = spec[k]
        if c in " |\n\t":
            k += 1
        elif c in "@!<":
            _, k = read_group(spec, k + 1)
        elif c == ">":
            g, k = read_group(spec, k + 1)
            pending = "c" if "\\centering" in g else "r" if "\\raggedleft" in g else "l"
        elif c == "*":
            count, k = read_group(spec, k + 1)
            sub, k = read_group(spec, k)
            cols += parse_colspec(sub, line_cm) * int(count)
        elif c in "lcr":
            cols.append((pending or c, None))
            pending, k = None, k + 1
        elif c in "pmbLRC":
            w, k = read_group(spec, k + 1)
            cols.append((pending or {"R": "r", "C": "c"}.get(c, "l"), width_cm(w, line_cm)))
            pending = None
        elif c in "XY":
            cols.append((pending or "l", None))
            pending, k = None, k + 1
        else:
            raise LatexError(f"column type {c!r} in {spec!r}")
    return cols


def colgroup(cols: list[tuple[str, float | None]], flexible: bool, line_cm: float) -> str:
    """<colgroup> with the LaTeX column widths as percentages: of the line width when the table has flexible (X)
    columns, of the fixed widths' sum otherwise (a tabular is as wide as its columns)."""
    fixed = [w for _, w in cols if w]
    if not fixed:
        return ""
    total = line_cm if flexible else sum(fixed)
    return "<colgroup>" + "".join(f"<col style=\"width:{100 * w / total:.1f}%\">" if w else "<col>"
                                  for _, w in cols) + "</colgroup>"


def split_top(raw: str, row: bool) -> list[str]:
    """Split a tabular body at top-level \\\\ (rows, with an optional [skip]) or & (cells)."""
    parts, depth, start, k, in_math = [], 0, 0, 0, False
    while k < len(raw):
        c = raw[k]
        if c == "\\":
            if row and raw.startswith("\\\\", k) and depth == 0 and not in_math:
                parts.append(raw[start:k])
                k += 2
                m = re.match(r"\s*\[[^\]]*\]", raw[k:])
                if m:
                    k += m.end()
                start = k
                continue
            k += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
        elif c == "$":
            in_math = not in_math
        elif c == "&" and not row and depth == 0 and not in_math:
            parts.append(raw[start:k])
            start = k + 1
        k += 1
    parts.append(raw[start:])
    return parts


RULE = re.compile(r"\s*\\(toprule|midrule|bottomrule|hline|cmidrule|addlinespace|morecmidrules)\b"
                  r"(?:\s*\[[^\]]*\])?(?:\s*\(([^)]*)\))?(?:\s*\{(\d+)-(\d+)\})?")


def table_rows(raw: str) -> list[dict]:
    rows: list[dict] = []
    pending: list[str] = []
    for chunk in split_top(raw, row=True):
        while True:
            m = RULE.match(chunk)
            if not m:
                break
            if m.group(1) == "cmidrule" and rows:
                rows[-1]["cmid_below"].append((int(m.group(3)), int(m.group(4))))
            else:
                pending.append(m.group(1))
            chunk = chunk[m.end():]
        if not chunk.strip():
            continue
        cells = []
        for cell in split_top(chunk, row=False):
            cell = cell.strip()
            m = re.match(r"\\multicolumn\s*", cell)
            if m:
                p = Parser(cell)
                p.i = m.end()
                span, spec = p.mandatory(), p.mandatory()
                content = p.mandatory()
                aligns = parse_colspec(spec)
                cells.append((int(span), aligns[0][0] if aligns else "l", content))
            else:
                cells.append((1, None, cell))
        rows.append({"cells": cells, "rules": pending, "cmid_below": []})
        pending = []
    return rows


# ------------------------------------------------------------------------------------------------------ page
# Charter (macOS) draws U+00A0 twice as wide as a space; the "Paper Serif" alias takes every other character from
# Charter and leaves the no-break space to the next font of the stack.
CSS = r"""
@font-face{font-family:"Paper Serif";src:local("Charter Roman"),local("Charter-Roman");font-weight:400;
font-style:normal;unicode-range:U+0-9F,U+A1-10FFFF}
@font-face{font-family:"Paper Serif";src:local("Charter Italic"),local("Charter-Italic");font-weight:400;
font-style:italic;unicode-range:U+0-9F,U+A1-10FFFF}
@font-face{font-family:"Paper Serif";src:local("Charter Bold"),local("Charter-Bold");font-weight:700;
font-style:normal;unicode-range:U+0-9F,U+A1-10FFFF}
@font-face{font-family:"Paper Serif";src:local("Charter Bold Italic"),local("Charter-BoldItalic");font-weight:700;
font-style:italic;unicode-range:U+0-9F,U+A1-10FFFF}
:root{--text:#1d1d1f;--muted:#5d5d63;--rule:#cfcfd4;--strong-rule:#1d1d1f;--link:#0a5292;--hl:#fff3cf;--bg:#fff;
--serif:"Paper Serif","Iowan Old Style","Bitstream Charter","Sitka Text",Cambria,Georgia,serif;
--sans:-apple-system,BlinkMacSystemFont,"Segoe UI","Helvetica Neue",Arial,sans-serif;color-scheme:light}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%;text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--text);font-family:var(--serif);font-size:1.0625rem;line-height:1.6;
font-kerning:normal;text-rendering:optimizeLegibility}
a{color:var(--link);text-decoration:none}
a:hover{text-decoration:underline}
.topline{border-bottom:1px solid var(--rule);font:0.8125rem/1.4 var(--sans);color:var(--muted)}
.topline .inner{max-width:46em;margin:0 auto;padding:.7em 16px;display:flex;flex-wrap:wrap;gap:.3em 1.2em;
justify-content:space-between}
.topline .journal{font-weight:600;color:var(--text)}
.topline nav .sep{margin:0 .55em;color:var(--rule)}
main{max-width:46em;margin:0 auto;padding:0 16px 3em}
header.front{padding:1.6em 0 .6em}
.kicker{font:600 .8rem/1.3 var(--sans);letter-spacing:.06em;text-transform:uppercase;color:var(--muted);margin:0 0 .5em}
h1.title{font-size:1.85em;line-height:1.22;font-weight:normal;margin:0 0 .7em;letter-spacing:-.005em}
.authors{font-size:1.06em;margin:0 0 .5em}
.authors sup,.affiliations sup,.corresp sup{font-size:.7em;margin-left:.08em}
.affiliations{list-style:none;padding:0;margin:0 0 .6em;font-size:.86em;font-style:italic;color:var(--muted);
line-height:1.45}
.affiliations sup{margin-right:.25em;font-style:normal}
.corresp{font-size:.86em;color:var(--muted);margin:0}
h2{font-size:1.24em;line-height:1.3;margin:2.1em 0 .65em;font-weight:bold}
h3{font-size:1.06em;line-height:1.35;margin:1.7em 0 .5em;font-weight:normal;font-style:italic}
h4{font-size:1em;margin:1.4em 0 .4em;font-weight:normal;font-style:italic}
.secnum{margin-right:.25em}
p{margin:0 0 .85em}
section.highlights,section.abstract{border-top:1px solid var(--strong-rule);margin-top:1.4em}
section.highlights h2,section.abstract h2{font-size:1.06em;margin:.8em 0 .5em}
section.highlights ul{margin:0 0 1em;padding-left:1.3em}
section.highlights li{margin:.2em 0}
section.abstract{border-bottom:1px solid var(--strong-rule);padding-bottom:.4em;margin-bottom:1.6em}
.keywords{font-size:.94em}
.kw-label{font-style:italic;margin-right:.25em}
.math.display{flex:1 1 auto;min-width:0;overflow-x:auto;overflow-y:hidden}
.equation{display:flex;align-items:center;gap:1em;margin:.9em 0}
.eqno{flex:none}
figure{margin:2em 0}
figure.figure img{display:block;width:100%;height:auto;margin:0 auto}
.img-link{display:block}
figcaption{font-size:.87em;line-height:1.5;margin-top:.7em}
figcaption .label{font-weight:bold}
figure.table figcaption{margin:0 0 .6em}
.table-wrap{overflow-x:auto;-webkit-overflow-scrolling:touch}
table{border-collapse:collapse;width:100%;font-size:.84em;line-height:1.4;font-variant-numeric:lining-nums tabular-nums;
border-top:1.5px solid var(--strong-rule);border-bottom:1.5px solid var(--strong-rule)}
th,td{padding:.32em .55em;vertical-align:top;text-align:left}
th{font-weight:normal;vertical-align:bottom}
thead tr:last-child th{border-bottom:1px solid var(--strong-rule)}
th.cmid{border-bottom:.5px solid var(--strong-rule)}
td.c,th.c{text-align:center}
td.r,th.r{text-align:right}
tr.midrule td{border-top:1px solid var(--strong-rule)}
tr.group td{font-style:normal;padding-top:.55em}
.tablenotes{font-size:.8em;line-height:1.45;margin-top:.55em;color:#333}
figure.cols5 table,figure.cols6 table{min-width:34em}
figure.cols7 table,figure.cols8 table{min-width:44em}
@media (min-width:62em){figure.table.wide{width:min(64em,calc(100vw - 32px));margin-left:50%;
transform:translateX(-50%)}}
.references h2{margin-top:2.4em}
.reflist{list-style:none;padding:0;margin:0;font-size:.9em;line-height:1.5}
.reflist li{padding-left:1.4em;text-indent:-1.4em;margin:0 0 .5em}
.reflist a{overflow-wrap:anywhere}
a.url{overflow-wrap:anywhere}
:target{scroll-margin-top:.8em}
li:target,figure:target>figcaption,.equation:target,h2:target,h3:target,h4:target{background:var(--hl);
box-shadow:0 0 0 .35em var(--hl);border-radius:2px}
.sc{font-variant:small-caps}
footer.colophon{max-width:46em;margin:0 auto;padding:1.2em 16px 3em;border-top:1px solid var(--rule);
font:.78rem/1.5 var(--sans);color:var(--muted)}
mjx-container[jax="CHTML"]{line-height:0}
@media (max-width:600px){body{font-size:1rem}h1.title{font-size:1.45em}.equation{gap:.5em}}
@media print{.topline,footer.colophon{display:none}a{color:inherit}main{max-width:none}}
"""

MATHJAX_CONFIG = r"""window.MathJax={tex:{inlineMath:[['\\(','\\)']],displayMath:[['\\[','\\]']],
macros:{symbf:['\\mathbf{#1}',1]}},chtml:{matchFontHeight:true},options:{enableMenu:false}};"""


def page_html(title: str, description: str, topnav: str, body: str, colophon: str, journal: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<meta name="description" content="{attr(description)}">
<style>{CSS}</style>
<script>{MATHJAX_CONFIG}</script>
<script defer src="{MATHJAX_URL}"></script>
</head>
<body>
<div class="topline"><div class="inner"><span class="journal">{esc(journal)}</span><nav>{topnav}</nav></div></div>
<main>
{body}
</main>
<footer class="colophon">{colophon}</footer>
</body>
</html>
"""


# ------------------------------------------------------------------------------------------------------ driver
def compile_intermediates(repo: Path, work: Path) -> Path:
    """Copy the repository's rendered sources into work/build and compile them there, keeping .aux and .bbl."""
    build = work / "build"
    build.mkdir(parents=True, exist_ok=True)
    figures = work / "figures"
    if not figures.exists():
        figures.symlink_to(repo / "paper" / "figures")
    bib_src, bib_dst = repo / "paper" / "references.bib", build / "references.bib"
    bib_same = bib_dst.exists() and bib_dst.read_bytes() == bib_src.read_bytes()
    shutil.copy2(bib_src, bib_dst)
    for doc in DOCS:
        src = repo / "paper" / "build" / f"{doc}.tex"
        dst = build / f"{doc}.tex"
        kept = [build / f"{doc}.{ext}" for ext in ("aux", "bbl", "pdf")]
        fresh = all(k.exists() and k.stat().st_mtime >= src.stat().st_mtime for k in kept)
        if fresh and bib_same and dst.exists() and dst.read_bytes() == src.read_bytes():
            continue
        shutil.copy2(src, dst)
        t = time.time()
        r = subprocess.run(["tectonic", "-X", "compile", "--only-cached", "--keep-intermediates", "--keep-logs",
                            dst.name], cwd=build, capture_output=True, text=True)
        if r.returncode:
            raise SystemExit(f"tectonic failed for {doc}:\n{(r.stdout + r.stderr)[-3000:]}")
        print(f"compiled {doc}.tex with tectonic in {time.time() - t:.1f} s (kept .aux, .bbl)")
    return build


def label_links(aux_labels: dict, page: str) -> dict:
    """(kind, number) -> URL of every section, figure, table and equation label of a document."""
    links = {}
    for key, number in aux_labels.items():
        kind = key.split(":", 1)[0]
        if kind in ("sec", "fig", "tab", "eq", "app") and number:
            links.setdefault(("sec" if kind == "app" else kind, number), f"{page}#{html_id(key)}")
    return links


def build(build_dir: Path, figures: Path, out: Path, dpi: int, docs: list[str], cross_links: bool) -> int:
    t0 = time.time()
    out.mkdir(parents=True, exist_ok=True)
    images = out / "images"
    images.mkdir(exist_ok=True)
    aux = {doc: read_aux(build_dir / f"{doc}.aux") for doc in DOCS if (build_dir / f"{doc}.aux").exists()}
    jobs: list[ImageJob] = []
    pages = {}
    m = re.search(r"\\journal\{([^}]*)\}", (build_dir / "manuscript.tex").read_text(encoding="utf-8"))
    journal_main = m.group(1) if m else ""
    for doc in docs:
        tex_path = build_dir / f"{doc}.tex"
        tex = tex_path.read_text(encoding="utf-8")
        if doc not in aux:
            raise SystemExit(f"{build_dir / (doc + '.aux')} not found: compile with tectonic --keep-intermediates")
        for ext in ("aux", "bbl"):
            if (build_dir / f"{doc}.{ext}").stat().st_mtime < tex_path.stat().st_mtime:
                print(f"  warning: {doc}.{ext} is older than {doc}.tex: numbers and references may be stale")
        bib = read_bbl(build_dir / f"{doc}.bbl")
        other = "supplementary" if doc == "manuscript" else "manuscript"
        links = label_links(aux[other][0], DOCS[other]) if (cross_links and other in aux) else {}
        if doc == "manuscript":  # only S-numbered items point to the supplement
            links = {k: v for k, v in links.items() if k[1].startswith("S")}
        else:
            links = {k: v for k, v in links.items() if not k[1].startswith(("S", "Appendix"))}
        r = Renderer(doc, tex, aux[doc], bib, figures, images, links)
        body = r.render()
        jobs += [j for j in r.jobs if all(j.png != k.png for k in jobs)]
        pages[doc] = (r, body)
    with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 2)) as pool:
        list(pool.map(lambda j: rasterize(j, dpi), jobs))
    sizes = {f"images/{j.png.name}": png_size(j.png) for j in jobs}
    for doc, (r, body) in pages.items():
        body = re.sub(r'<img src="(images/[^"]+)" width="0" height="0"',
                      lambda m: '<img src="{}" width="{}" height="{}"'.format(m.group(1), *sizes[m.group(1)]), body)
        if doc == "manuscript":
            nav = ("<a href=\"manuscript.pdf\">PDF</a><span class=\"sep\">|</span>Supplementary material: "
                   "<a href=\"supplementary.html\">HTML</a> · <a href=\"supplementary.pdf\">PDF</a>")
            abstract = re.search(r"<section class=\"abstract\"[^>]*><h2>Abstract</h2><p>(.*?)</p>", body, re.S)
            description = strip_tags(abstract.group(1))[:300] if abstract else ""
        else:
            nav = ("<a href=\"supplementary.pdf\">PDF</a><span class=\"sep\">|</span>Main text: "
                   "<a href=\"index.html\">HTML</a> · <a href=\"manuscript.pdf\">PDF</a>")
            description = "Supplementary material of the manuscript."
        colophon = ("HTML rendering generated from the LaTeX source of the "
                    + ("manuscript" if doc == "manuscript" else "supplementary material")
                    + "; equations typeset by MathJax. The PDF is the reference version.")
        journal = r.journal or journal_main
        kicker = "" if doc == "manuscript" else "Supplementary material · "
        html_text = page_html(r.title_text or doc, description, nav, body, colophon,
                              kicker + (f"Preprint submitted to {journal}" if journal else "Preprint"))
        (out / DOCS[doc]).write_text(html_text, encoding="utf-8")
        pdf = build_dir / f"{doc}.pdf"
        if pdf.exists():
            shutil.copy2(pdf, out / f"{doc}.pdf")
        checked = sum(1 for k in r.labels if k in r.aux_labels)
        print(f"wrote {out / DOCS[doc]} ({len(html_text) / 1024:.0f} KiB): {len(r.labels)} labels "
              f"({checked} checked against the .aux), {len(r.bib)} references, {r.n_cites} citations, "
              f"{len(r.jobs)} images")
        for w in r.warnings:
            print(f"  warning: {w}")
    if set(docs) == set(DOCS):  # images of figures no longer in either page
        used = {j.png.name for j in jobs}
        for stale in images.glob("*.png"):
            if stale.name not in used:
                stale.unlink()
    n_warn = sum(len(r.warnings) for r, _ in pages.values())
    print(f"{len(jobs)} images at {dpi} dpi; {n_warn} warnings; done in {time.time() - t0:.1f} s")
    return n_warn



# ------------------------------------------------------------------------------------------------------ link check
class _Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.refs, self.ids, self.dup = [], set(), []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if a.get("id"):
            if a["id"] in self.ids:
                self.dup.append(a["id"])
            self.ids.add(a["id"])
        for key in ("href", "src"):
            if a.get(key):
                self.refs.append(a[key])


def check_links(root: Path) -> list[str]:
    """Broken local links in every HTML file under ``root``: missing files, missing #anchors (also within a page:
    cross-references, the reference list) and duplicate ids. External links (http, https, mailto) are not checked
    here; scripts/check_live_site.py checks the published pages."""
    root = root.resolve()
    pages, problems = {}, []
    for f in root.rglob("*.html"):
        p = _Links()
        p.feed(f.read_text(encoding="utf-8"))
        pages[f.resolve()] = p
    for f, p in pages.items():
        problems += [f"{f.relative_to(root)}: duplicate id {i}" for i in p.dup]
        for ref in p.refs:
            if re.match(r"^(https?:|mailto:|data:)", ref):
                continue
            target, _, anchor = ref.partition("#")
            path = (f.parent / target).resolve() if target else f
            if not path.exists():
                problems.append(f"{f.relative_to(root)}: missing {ref}")
            elif anchor and path.suffix == ".html" and anchor not in pages.get(path, _Links()).ids:
                problems.append(f"{f.relative_to(root)}: missing anchor {ref}")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", type=Path, help="repository to read paper/build/*.tex, references.bib and figures from")
    ap.add_argument("--work", type=Path, default=Path("texbuild"), help="scratch folder for the tectonic run (--repo)")
    ap.add_argument("--build", type=Path, help="folder holding <doc>.tex, .aux, .bbl and .pdf (instead of --repo)")
    ap.add_argument("--figures", type=Path, help="folder of Figure_*.pdf (default: <repo>/paper/figures)")
    ap.add_argument("--out", type=Path, default=Path("html"))
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--only", choices=list(DOCS), help="build one page")
    ap.add_argument("--no-cross-links", action="store_true", help="do not link 'Fig. S6' etc. to the other page")
    ap.add_argument("--strict", action="store_true", help="exit with an error when the conversion warns")
    args = ap.parse_args(argv)
    if args.repo:
        build_dir = compile_intermediates(args.repo.resolve(), args.work.resolve())
        figures = args.figures or args.repo / "paper" / "figures"
    elif args.build:
        build_dir = args.build
        figures = args.figures or build_dir.parent / "figures"
    else:
        ap.error("give --repo or --build")
    docs = [args.only] if args.only else list(DOCS)
    n_warn = build(build_dir, figures, args.out, args.dpi, docs, not args.no_cross_links)
    broken = check_links(args.out)
    for b in broken:
        print(f"  broken link: {b}")
    if args.strict and (n_warn or broken):
        raise SystemExit(f"{n_warn} warnings and {len(broken)} broken links (--strict)")


if __name__ == "__main__":
    main()
