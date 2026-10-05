"""Helpers for the static site (G5): a small Markdown subset to HTML, figures, tables, the report's
number facts, file manifests and a link checker. The page content lives in scripts/build_site.py.

The Markdown subset is what the project documents use: ATX headings (an explicit id may end the
line: ``## 2. Methods {#methods}``), paragraphs, bullet and numbered lists with indented
continuation lines (one nesting level), pipe tables, fenced code, inline code, bold, italic and
links. For the report also: figure blocks (``::: figure {#fig-id}``, then ``![alt](results/....png)``
lines and the caption, closed by ``:::``), table captions (a line ``Table: caption {#tab-id}``
before a pipe table), both numbered in order of appearance, cross-references to them (``[@fig-id]``,
``[@tab-id]``), subscripts (``X~child~``) and superscripts (``cm^2^``). Raw HTML is always escaped.

Two documents can refer to each other's figures and tables (the report and its supplementary text):
``md_blocks`` converts a document and returns its labels (with a number prefix, e.g. 'S' for
'Figure S1'), ``resolve_xrefs`` then links each ``[@id]`` to its own document ('#id') or to the
other ('page.html#id'); ``md_to_html`` does both for a single document.
"""
from __future__ import annotations

import hashlib
import html
import re
import shutil
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath

import jinja2

_LIST = re.compile(r"^(\s*)([*-]|\d+\.)\s+(.*)$")
_ID = re.compile(r"^[A-Za-z][\w:.-]*$")
_END_ID = re.compile(r"\s*\{#([^{}]*)\}\s*$")
_FIGURE = re.compile(r"^:::\s*figure\s*(?:\{#([^{}]*)\})?\s*$")
_IMAGE = re.compile(r"^\s*!\[([^\]]*)\]\(([^)\s]+)\)\s*$")
# sub- and superscripts: no space inside and no letter or digit right after (b_k^2/s_k^2 stays as written)
_SUB = re.compile(r"~([^~\s]+)~(?!\w)")
_SUP = re.compile(r"\^([^^\s]+)\^(?!\w)")
_XREF = re.compile(r"\[@([^\]\s]+)\]")


def slug(text: str) -> str:
    s = re.sub(r"<[^>]+>", "", text).lower()
    s = re.sub(r"[^a-z0-9]+", "-", s).strip("-")
    return s or "section"


def _outside_tags(h: str, fn) -> str:
    """Apply ``fn`` to the text between the tags of ``h`` (escaped text: every '<' starts a tag)."""
    parts = re.split(r"(<[^>]*>)", h)
    return "".join(p if k % 2 else fn(p) for k, p in enumerate(parts))


def inline(text: str) -> str:
    """Escape HTML, then code spans, bold, italic, links, sub- and superscripts (never inside code spans)."""
    out = []
    for part in re.split(r"(`[^`]+`)", text):
        if len(part) > 1 and part.startswith("`") and part.endswith("`"):
            out.append(f"<code>{html.escape(part[1:-1], quote=False)}</code>")
            continue
        p = html.escape(part, quote=False)
        p = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", p)
        p = re.sub(r"(?<![\w*])\*(?![\s*])(.+?)(?<![\s*])\*(?![\w*])", r"<em>\1</em>", p)
        p = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)",
                   lambda m: f'<a href="{html.escape(html.unescape(m.group(2)))}">{m.group(1)}</a>', p)
        out.append(_outside_tags(p, lambda t: _SUP.sub(r"<sup>\1</sup>", _SUB.sub(r"<sub>\1</sub>", t))))
    return "".join(out)


def _valid_id(ident: str) -> str:
    if not _ID.match(ident):
        raise ValueError(f"invalid id {{#{ident}}}: a letter, then letters, digits, '-', '_', ':' or '.'")
    return ident


def _split_id(text: str) -> tuple[str, str | None]:
    """'Text {#id}' -> ('Text', 'id'); raises on an invalid id."""
    m = _END_ID.search(text)
    return (text[:m.start()], _valid_id(m.group(1))) if m else (text, None)


def _is_table(lines: list[str], i: int) -> bool:
    return lines[i].lstrip().startswith("|") and i + 1 < len(lines) and bool(re.match(r"^\s*\|?\s*:?-{2,}", lines[i + 1]))


def _table(rows: list[str], caption_html: str | None = None, ident: str | None = None) -> str:
    def cells(line):
        return [c.strip() for c in line.strip().strip("|").split("|")]

    head, body = cells(rows[0]), [cells(r) for r in rows[2:]]
    h = "".join(f"<th>{inline(c)}</th>" for c in head)
    b = "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in body)
    attr = f' id="{ident}"' if ident else ""
    cap = f"<caption>{caption_html}</caption>" if caption_html else ""
    return f'<div class="table-wrap"><table{attr}>{cap}<thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def md_to_html(md: str, heading_offset: int = 0, image=None, prefix: str = "", refs: dict | None = None) -> str:
    """Convert the project's Markdown subset to HTML (``heading_offset`` demotes headings).

    ``image`` maps an image path of a figure block to its published src (e.g. publish_png); without
    it the path is used as written. ``prefix`` goes before the figure and table numbers ('S':
    'Figure S1'); ``refs`` are the figures and tables of another document (see resolve_xrefs).
    Raises ValueError on a malformed block, a duplicate or invalid id, or a cross-reference to an
    unknown id.
    """
    return resolve_xrefs(*md_blocks(md, heading_offset, image, prefix), refs)


def md_blocks(md: str, heading_offset: int = 0, image=None, prefix: str = "") -> tuple[str, dict]:
    """The conversion of md_to_html without the cross-references: returns the HTML, in which every
    ``[@id]`` is left as written, and the document's figure and table labels, id -> 'Figure 1' (with
    ``prefix``: 'Figure S1'), in order of appearance."""
    lines = md.splitlines()
    out, para, i = [], [], 0
    ids, labels, count = set(), {}, {"Figure": 0, "Table": 0}

    def flush():
        if para:
            out.append(f"<p>{inline(' '.join(x.strip() for x in para))}</p>")
            para.clear()

    def claim(ident, auto=None):
        """Reserve an explicit id (a duplicate raises) or a unique heading slug."""
        if ident is None:
            ident, k = auto, 2
            while ident in ids:
                ident, k = f"{auto}-{k}", k + 1
        elif ident in ids:
            raise ValueError(f"duplicate id {ident!r}")
        ids.add(ident)
        return ident

    def number(kind, ident):
        count[kind] += 1
        label = f"{kind} {prefix}{count[kind]}"
        if ident is not None:
            labels[claim(ident)] = label
        return label

    while i < len(lines):
        line = lines[i]
        if not line.strip():
            flush()
            i += 1
            continue
        if line.lstrip().startswith("<!--"):
            flush()
            while i < len(lines) and "-->" not in lines[i]:
                i += 1
            i += 1
            continue
        if line.startswith("```"):
            flush()
            j = i + 1
            while j < len(lines) and not lines[j].startswith("```"):
                j += 1
            out.append("<pre><code>" + html.escape("\n".join(lines[i + 1:j]), quote=False) + "</code></pre>")
            i = j + 1
            continue
        if line.startswith(":::"):
            flush()
            m = _FIGURE.match(line)
            if not m:
                raise ValueError(f"line {i + 1}: unsupported block {line!r} (only '::: figure {{#fig-id}}' ... ':::')")
            ident = _valid_id(m.group(1)) if m.group(1) is not None else None
            j = i + 1
            while j < len(lines) and not lines[j].lstrip().startswith(":::"):
                j += 1
            if j == len(lines) or lines[j].strip() != ":::":
                raise ValueError(f"line {i + 1}: figure block not closed by ':::'")
            images, caption = [], []
            for x in lines[i + 1:j]:
                im = _IMAGE.match(x)
                if im:
                    images.append(im.groups())
                elif x.strip():
                    caption.append(x.strip())
            if not images:
                raise ValueError(f"line {i + 1}: figure block without an image")
            label = number("Figure", ident)
            imgs = []
            for alt, path in images:
                src = html.escape(image(path) if image else path)
                imgs.append(f'<a href="{src}"><img src="{src}" alt="{html.escape(alt)}" loading="lazy"></a>')
            cap = f"<strong>{label}.</strong>" + (f" {inline(' '.join(caption))}" if caption else "")
            attr = f' id="{ident}"' if ident else ""
            out.append(f'<figure{attr}><div class="figure-images">{"".join(imgs)}</div>'
                       f"<figcaption>{cap}</figcaption></figure>")
            i = j + 1
            continue
        if line.startswith("Table:"):
            flush()
            text, ident = _split_id(line[len("Table:"):].strip())
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j == len(lines) or not _is_table(lines, j):
                raise ValueError(f"line {i + 1}: a table caption must be followed by a pipe table")
            k = j
            while k < len(lines) and lines[k].lstrip().startswith("|"):
                k += 1
            label = number("Table", ident)
            out.append(_table(lines[j:k], f"<strong>{label}.</strong> {inline(text)}", ident))
            i = k
            continue
        m = re.match(r"^(#{1,6})\s+(.*)$", line)
        if m:
            flush()
            level = min(6, len(m.group(1)) + heading_offset)
            raw, ident = _split_id(m.group(2).strip())
            text = inline(raw)
            out.append(f'<h{level} id="{claim(ident, slug(text))}">{text}</h{level}>')
            i += 1
            continue
        if _is_table(lines, i):
            flush()
            j = i
            while j < len(lines) and lines[j].lstrip().startswith("|"):
                j += 1
            out.append(_table(lines[i:j]))
            i = j
            continue
        m = _LIST.match(line)
        if m and len(m.group(1)) == 0:
            flush()
            ordered = m.group(2)[0].isdigit()
            items = []  # (text, [nested texts])
            while i < len(lines):
                m = _LIST.match(lines[i])
                if m and len(m.group(1)) == 0 and m.group(2)[0].isdigit() == ordered:
                    items.append([m.group(3).strip(), []])
                elif m and len(m.group(1)) > 0 and items:
                    items[-1][1].append(m.group(3).strip())
                elif lines[i].startswith((" ", "\t")) and lines[i].strip() and items:
                    if items[-1][1]:
                        items[-1][1][-1] += " " + lines[i].strip()
                    else:
                        items[-1][0] += " " + lines[i].strip()
                else:
                    break
                i += 1
            tag = "ol" if ordered else "ul"
            lis = []
            for text, nested in items:
                sub = "<ul>" + "".join(f"<li>{inline(t)}</li>" for t in nested) + "</ul>" if nested else ""
                lis.append(f"<li>{inline(text)}{sub}</li>")
            out.append(f"<{tag}>{''.join(lis)}</{tag}>")
            continue
        para.append(line)
        i += 1
    flush()
    return "\n".join(out), labels


def resolve_xrefs(h: str, labels: dict, refs: dict | None = None) -> str:
    """Replace every ``[@id]`` of ``h`` outside code and tags (after the conversion: references may point
    forward) by a link: to '#id' for a figure or table of this document (``labels``, id -> label), to
    'page#id' for one of another document (``refs``, id -> (page, label)). An unknown id, or an id that
    both documents define, raises ValueError."""
    refs = refs or {}
    both = sorted(set(labels) & set(refs))
    if both:
        raise ValueError("figure or table id defined in both documents: " + ", ".join(both))

    def ref(m):
        ident = m.group(1)
        if ident in labels:
            return f'<a href="#{ident}">{labels[ident]}</a>'
        if ident in refs:
            page, label = refs[ident]
            return f'<a href="{html.escape(page)}#{ident}">{label}</a>'
        raise ValueError(f"cross-reference to an unknown figure or table: [@{ident}]")

    parts = re.split(r"(<code>[\s\S]*?</code>|<[^>]*>)", h)
    return "".join(p if k % 2 else _XREF.sub(ref, p) for k, p in enumerate(parts))


def toc(h: str) -> str:
    """Put a table of contents of the h2/h3 headings of ``h`` before its first section heading."""
    tree = []  # [id, text, [(id, text)]]
    for level, ident, text in re.findall(r'<h([23]) id="([^"]+)">(.*?)</h\1>', h):
        text = re.sub(r"</?a\b[^>]*>", "", text)
        if level == "2" or not tree:
            tree.append((ident, text, []))
        else:
            tree[-1][2].append((ident, text))
    if not tree:
        return h

    def li(ident, text, sub=()):
        nested = "<ul>" + "".join(li(*s) for s in sub) + "</ul>" if sub else ""
        return f'<li><a href="#{ident}">{text}</a>{nested}</li>'

    nav = ('<nav class="toc" aria-label="Contents"><details open><summary>Contents</summary><ul>'
           + "".join(li(*t) for t in tree) + "</ul></details></nav>")
    k = re.search(r"<h[23] ", h).start()
    return h[:k] + nav + "\n" + h[k:]


class Facts:
    """``F`` in the report: ``F.name`` or ``F["name"]`` is the fact's value as printed; the names read
    are recorded in order. A missing name is undefined, which StrictUndefined turns into an error."""

    def __init__(self, values: dict):
        self._values, self._used = dict(values), []

    def __getitem__(self, name):
        value = self._values[name]
        if name not in self._used:
            self._used.append(name)
        return value


def fill_facts(text: str, values: dict) -> tuple[str, list[str]]:
    """Render ``text`` through Jinja2 with ``F`` (fact name -> value as printed) and StrictUndefined: a
    missing fact or any other variable raises jinja2.UndefinedError (with its line). Returns the text
    and the fact names used. Jinja comments are ``{## ... ##}``, so that ``{#id}`` stays Markdown."""
    env = jinja2.Environment(undefined=jinja2.StrictUndefined, autoescape=False, keep_trailing_newline=True,
                             comment_start_string="{##", comment_end_string="##}")
    facts = Facts(values)
    try:
        return env.from_string(text).render(F=facts), facts._used
    except jinja2.UndefinedError as e:
        tb, line = e.__traceback__, None
        while tb:  # Jinja maps template frames to their source lines
            line = tb.tb_lineno if tb.tb_frame.f_code.co_filename == "<template>" else line
            tb = tb.tb_next
        msg = re.sub(r"^'[\w.]*Facts object' has no attribute ", "no fact ", e.message or "")
        raise jinja2.UndefinedError(f"line {line}: {msg}" if line else msg) from None


def publish_png(path: str, results: Path, out: Path) -> str:
    """Copy a committed result figure, written 'results/<dir>/<name>.png', to out/figures/<dir>/<name>.png
    and return that href. Anything else (another folder, '..', another format, a missing file) raises."""
    p = PurePosixPath(path)
    if len(p.parts) < 2 or p.parts[0] != "results" or ".." in p.parts or p.suffix.lower() != ".png":
        raise ValueError(f"figure {path!r}: only PNG files under results/ are published")
    rel = PurePosixPath(*p.parts[1:])
    src = results / rel
    if not src.is_file():
        raise FileNotFoundError(f"figure {path} not found")
    dst = out / "figures" / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return f"figures/{rel}"


def table(header: list[str], rows: list[list], caption: str | None = None, cls: str = "", html_cols: tuple = ()) -> str:
    """HTML table from formatted cell strings, escaped except in ``html_cols`` (trusted markup)."""
    h = "".join(f"<th>{html.escape(str(c))}</th>" for c in header)
    b = "".join("<tr>" + "".join(f"<td>{c if j in html_cols else html.escape(str(c))}</td>" for j, c in enumerate(r)) + "</tr>"
                for r in rows)
    cap = f"<caption>{html.escape(caption)}</caption>" if caption else ""
    return f'<div class="table-wrap"><table class="{cls}">{cap}<thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def figure(src: str, caption_html: str, alt: str) -> str:
    return (f'<figure><a href="{html.escape(src)}"><img src="{html.escape(src)}" alt="{html.escape(alt)}" loading="lazy"></a>'
            f"<figcaption>{caption_html}</figcaption></figure>")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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
    """Broken local links in every HTML file under ``root``: missing files, missing #anchors (also
    within a page: tables of contents, cross-references) and duplicate ids. External (http, https,
    mailto) links are listed separately by the caller if needed."""
    root = root.resolve()
    pages, problems = {}, []
    for f in root.rglob("*.html"):
        p = _Links()
        p.feed(f.read_text())
        pages[f.resolve()] = p
    for f, p in pages.items():
        problems += [f"{f.relative_to(root)}: duplicate id {i}" for i in p.dup]
        for ref in p.refs:
            if re.match(r"^(https?:|mailto:)", ref):
                continue
            target, _, anchor = ref.partition("#")
            path = (f.parent / target).resolve() if target else f
            if not path.exists():
                problems.append(f"{f.relative_to(root)}: missing {ref}")
            elif anchor and path.suffix == ".html" and anchor not in pages.get(path, _Links()).ids:
                problems.append(f"{f.relative_to(root)}: missing anchor {ref}")
    return problems
