"""Helpers for the local static report (G5): a small Markdown subset to HTML, figures, tables,
file manifests and a link checker. The page content lives in scripts/build_site.py.

The Markdown subset is what the project documents use: ATX headings, paragraphs, bullet and
numbered lists with indented continuation lines (one nesting level), pipe tables, fenced code,
inline code, bold, italic and links.
"""
from __future__ import annotations

import hashlib
import html
import re
from html.parser import HTMLParser
from pathlib import Path

_LIST = re.compile(r"^(\s*)([*-]|\d+\.)\s+(.*)$")


def slug(text: str) -> str:
    s = re.sub(r"<[^>]+>", "", text).lower()
    s = re.sub(r"[^a-z0-9]+", "-", s).strip("-")
    return s or "section"


def inline(text: str) -> str:
    """Escape HTML, then code spans, bold, italic and links (never inside code spans)."""
    out = []
    for part in re.split(r"(`[^`]+`)", text):
        if len(part) > 1 and part.startswith("`") and part.endswith("`"):
            out.append(f"<code>{html.escape(part[1:-1], quote=False)}</code>")
            continue
        p = html.escape(part, quote=False)
        p = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", p)
        p = re.sub(r"(?<![\w*])\*(?![\s*])(.+?)(?<![\s*])\*(?![\w*])", r"<em>\1</em>", p)
        p = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", lambda m: f'<a href="{html.escape(m.group(2))}">{m.group(1)}</a>', p)
        out.append(p)
    return "".join(out)


def _table(rows: list[str]) -> str:
    def cells(line):
        return [c.strip() for c in line.strip().strip("|").split("|")]

    head, body = cells(rows[0]), [cells(r) for r in rows[2:]]
    h = "".join(f"<th>{inline(c)}</th>" for c in head)
    b = "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in body)
    return f'<div class="table-wrap"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def md_to_html(md: str, heading_offset: int = 0) -> str:
    """Convert the project's Markdown subset to HTML (``heading_offset`` demotes headings)."""
    lines = md.splitlines()
    out, para, i = [], [], 0

    def flush():
        if para:
            out.append(f"<p>{inline(' '.join(x.strip() for x in para))}</p>")
            para.clear()

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
        m = re.match(r"^(#{1,6})\s+(.*)$", line)
        if m:
            flush()
            level = min(6, len(m.group(1)) + heading_offset)
            text = inline(m.group(2).strip())
            out.append(f'<h{level} id="{slug(text)}">{text}</h{level}>')
            i += 1
            continue
        if line.lstrip().startswith("|") and i + 1 < len(lines) and re.match(r"^\s*\|?\s*:?-{2,}", lines[i + 1]):
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
    return "\n".join(out)


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
        self.refs, self.ids = [], set()

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if a.get("id"):
            self.ids.add(a["id"])
        for key in ("href", "src"):
            if a.get(key):
                self.refs.append(a[key])


def check_links(root: Path) -> list[str]:
    """Broken local links in every HTML file under ``root``: missing files or missing #anchors.
    External (http, https, mailto) links are listed separately by the caller if needed."""
    root = root.resolve()
    pages, problems = {}, []
    for f in root.rglob("*.html"):
        p = _Links()
        p.feed(f.read_text())
        pages[f.resolve()] = p
    for f, p in pages.items():
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
