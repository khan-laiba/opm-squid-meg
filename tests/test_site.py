import html
import importlib.util
import json
import re
import tempfile
import unittest
from pathlib import Path

import jinja2

from opmsquid import sitebuild as sb

ROOT = Path(__file__).resolve().parents[1]


def builder():
    """scripts/build_site.py as a module."""
    spec = importlib.util.spec_from_file_location("build_site", ROOT / "scripts" / "build_site.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestMarkdown(unittest.TestCase):
    def test_blocks_and_inline(self):
        md = ("# Title\n\nA *new* **bold** `x<y` [link](page.html#a).\n\n* item one\n  continued\n* item two\n  - nested\n\n"
              "1. first\n2. second\n\n| a | b |\n|---|---|\n| 1 | <script> |\n\n```bash\necho <hi>\n```\n")
        h = sb.md_to_html(md, heading_offset=1)
        self.assertIn('<h2 id="title">Title</h2>', h)
        self.assertIn("<em>new</em>", h)
        self.assertIn("<strong>bold</strong>", h)
        self.assertIn("<code>x&lt;y</code>", h)
        self.assertIn('<a href="page.html#a">link</a>', h)
        self.assertIn("<li>item one continued</li>", h)
        self.assertIn("<ul><li>nested</li></ul>", h)
        self.assertIn("<ol><li>first</li><li>second</li></ol>", h)
        self.assertIn("<td>&lt;script&gt;</td>", h)  # never raw HTML from documents
        self.assertIn("echo &lt;hi&gt;", h)

    def test_link_checker_finds_missing_files_and_anchors(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "a.html").write_text('<h2 id="x">x</h2><a href="b.html">b</a><a href="a.html#x">ok</a><a href="a.html#y">bad</a>'
                                      '<img src="f.png"><a href="https://example.org">ext</a>')
            problems = sb.check_links(d)
        self.assertEqual(sorted(problems), ["a.html: missing anchor a.html#y", "a.html: missing b.html", "a.html: missing f.png"])

    def test_link_checker_finds_anchors_within_a_page_and_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "figures").mkdir()
            (d / "figures" / "x.png").write_bytes(b"")
            (d / "r.html").write_text('<a href="#fig-a">Figure 1</a> <a href="#tab-b">Table 1</a><figure id="fig-a">'
                                      '<img src="figures/x.png"></figure><h2 id="s">s</h2><h3 id="s">t</h3>')
            problems = sb.check_links(d)
        self.assertEqual(sorted(problems), ["r.html: duplicate id s", "r.html: missing anchor #tab-b"])


class TestReportMarkdown(unittest.TestCase):
    """The report's additions to the Markdown subset."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.res, self.out = d / "results", d / "site"
        for rel in ("g9/a.png", "g9/b.png", "report/Figure_R1_x.png"):
            (self.res / rel).parent.mkdir(parents=True, exist_ok=True)
            (self.res / rel).write_bytes(b"\x89PNG\r\n\x1a\n" + rel.encode())
        (self.res / "g9" / "c.svg").write_text("<svg/>")

    def tearDown(self):
        self.tmp.cleanup()

    def render(self, md, heading_offset=0):
        return sb.md_to_html(md, heading_offset, image=lambda p: sb.publish_png(p, self.res, self.out))

    def test_figure_blocks(self):
        md = ("See [@fig-two] and [@fig-one].\n\n"
              "::: figure {#fig-one}\n![First *alt*](results/g9/a.png)\n![Second](results/report/Figure_R1_x.png)\n"
              "**Bold title.** Caption with `code` and\na second line.\n:::\n\n"
              "::: figure {#fig-two}\n![Only](results/g9/b.png)\n:::\n")
        h = self.render(md)
        self.assertIn('<figure id="fig-one"><div class="figure-images"><a href="figures/g9/a.png"><img src="figures/g9/a.png" '
                      'alt="First *alt*" loading="lazy"></a><a href="figures/report/Figure_R1_x.png"><img '
                      'src="figures/report/Figure_R1_x.png" alt="Second" loading="lazy"></a></div>', h)
        self.assertIn("<figcaption><strong>Figure 1.</strong> <strong>Bold title.</strong> Caption with <code>code</code> and a "
                      "second line.</figcaption></figure>", h)
        self.assertIn('<figure id="fig-two"><div class="figure-images"><a href="figures/g9/b.png">', h)
        self.assertIn("<figcaption><strong>Figure 2.</strong></figcaption>", h)
        self.assertIn('See <a href="#fig-two">Figure 2</a> and <a href="#fig-one">Figure 1</a>.', h)  # forward references
        for rel in ("g9/a.png", "g9/b.png", "report/Figure_R1_x.png"):  # copied to figures/<path relative to results/>
            self.assertEqual((self.out / "figures" / rel).read_bytes(), (self.res / rel).read_bytes())
        self.assertEqual(sb.publish_png("results/g9/a.png", self.res, self.out), "figures/g9/a.png")

    def test_figure_block_errors(self):
        img = "![x](results/g9/a.png)\n"
        bad = {"::: figure {#fig-a}\n" + img: ValueError,  # not closed
               "::: figure {#fig-a}\n" + img + "::: figure {#fig-b}\n![y](results/g9/b.png)\n:::\n": ValueError,
               "::: figure {#fig-a}\nA caption only.\n:::\n": ValueError,  # no image
               "::: figure {#fig-a}\n![x](docs/a.png)\n:::\n": ValueError,  # only results/
               "::: figure {#fig-a}\n![x](results/../a.png)\n:::\n": ValueError,
               "::: figure {#fig-a}\n![x](results/g9/c.svg)\n:::\n": ValueError,  # only PNG
               "::: figure {#fig-a}\n![x](results/g9/missing.png)\n:::\n": FileNotFoundError,
               "::: figure {#1a}\n" + img + ":::\n": ValueError,  # invalid id
               "::: note\ntext\n:::\n": ValueError}  # unsupported block
        for md, err in bad.items():
            with self.subTest(md=md), self.assertRaises(err):
                self.render(md)

    def test_table_captions(self):
        md = ("Table: First *caption* {#tab-a}\n\n| a | b |\n|---|---|\n| 1 | 2 |\n\n"
              "| plain | table |\n|---|---|\n| x | y |\n\n"
              "Table: Second caption\n| c |\n|---|\n| 3 |\n\nSee [@tab-a].\n")
        h = self.render(md)
        self.assertIn('<table id="tab-a"><caption><strong>Table 1.</strong> First <em>caption</em></caption>'
                      '<thead><tr><th>a</th>', h)
        self.assertIn("<table><thead><tr><th>plain</th>", h)  # a table without a caption keeps no number
        self.assertIn("<table><caption><strong>Table 2.</strong> Second caption</caption><thead>", h)
        self.assertIn('See <a href="#tab-a">Table 1</a>.', h)
        with self.assertRaises(ValueError):
            self.render("Table: a caption without a table\n\nA paragraph.\n")

    def test_cross_references(self):
        md = ("::: figure {#fig-a}\n![x](results/g9/a.png)\nCaption, see [@tab-a].\n:::\n\n"
              "Table: T {#tab-a}\n| a |\n|---|\n| [@fig-a] |\n\n"
              "[@fig-a] and `[@fig-a]` in code; *see [@tab-a]*.\n\n* item [@tab-a]\n")
        h = self.render(md)
        self.assertEqual(h.count('<a href="#fig-a">Figure 1</a>'), 2)
        self.assertEqual(h.count('<a href="#tab-a">Table 1</a>'), 3)
        self.assertIn("<code>[@fig-a]</code>", h)  # never inside code
        for md in ("See [@fig-missing].", "::: figure {#fig-a}\n![x](results/g9/a.png)\n:::\n\nSee [@tab-a].",
                   "::: figure {#fig-a}\n![x](results/g9/a.png)\n:::\n\n::: figure {#fig-a}\n![x](results/g9/b.png)\n:::\n",
                   "## Results {#fig-a}\n\n::: figure {#fig-a}\n![x](results/g9/a.png)\n:::\n"):  # unknown or duplicate ids
            with self.subTest(md=md), self.assertRaises(ValueError):
                self.render(md)

    def test_heading_ids(self):
        h = sb.md_to_html("# Title\n\n## 2. Methods {#methods}\n\n### Notes\n\n### Notes\n")
        self.assertIn('<h1 id="title">Title</h1>', h)
        self.assertIn('<h2 id="methods">2. Methods</h2>', h)
        self.assertIn('<h3 id="notes">Notes</h3>', h)
        self.assertIn('<h3 id="notes-2">Notes</h3>', h)  # a repeated heading gets a unique id
        self.assertIn('<h3 id="methods">2. Methods</h3>', sb.md_to_html("## 2. Methods {#methods}", heading_offset=1))
        for md in ("## A {#1bad}", "## A {#a b}", "## A {#x}\n\n## B {#x}"):
            with self.subTest(md=md), self.assertRaises(ValueError):
                sb.md_to_html(md)

    def test_sub_and_superscripts(self):
        h = sb.inline("D~child~ - D~adult~, 20 cm^2^, `X~a~ cm^2^`, [x~1~](https://example.org/~u/a~b?p=1&q=2), "
                      "b_k^2/s_k^2 and ~5 mm")
        self.assertIn("D<sub>child</sub> - D<sub>adult</sub>", h)
        self.assertIn("20 cm<sup>2</sup>", h)
        self.assertIn("<code>X~a~ cm^2^</code>", h)  # not inside code spans
        self.assertIn('<a href="https://example.org/~u/a~b?p=1&amp;q=2">x<sub>1</sub></a>', h)  # nor in link targets
        self.assertIn("b_k^2/s_k^2 and ~5 mm", h)  # formulas and approximations stay as written

    def test_table_of_contents(self):
        h = sb.toc(sb.md_to_html("Abstract.\n\n## 1. One {#one}\n\n### 1.1 Sub\n\n## 2. Two\n"))
        self.assertTrue(h.startswith('<p>Abstract.</p>\n<nav class="toc" aria-label="Contents">'))
        self.assertIn('<ul><li><a href="#one">1. One</a><ul><li><a href="#1-1-sub">1.1 Sub</a></li></ul></li>'
                      '<li><a href="#2-two">2. Two</a></li></ul>', h)
        self.assertEqual(sb.toc("<p>No sections.</p>"), "<p>No sections.</p>")


class TestFacts(unittest.TestCase):
    def test_fill_facts(self):
        text, used = sb.fill_facts("{{ F.b }} and {{ F['a'] }}, {{ F.b }}, {{ F.items }}.{## comment ##}\n"
                                   "## 2. Methods {#methods}\n",
                                   {"a": "1.14x", "b": "+0.44 dB", "items": "3", "c": "unused"})
        self.assertEqual(text, "+0.44 dB and 1.14x, +0.44 dB, 3.\n## 2. Methods {#methods}\n")  # {#id} is not a Jinja comment
        self.assertEqual(used, ["b", "a", "items"])

    def test_missing_fact_or_variable_raises(self):  # StrictUndefined
        for text in ("{{ F.missing }}", "{{ F['missing'] }}", "{{ other }}", "{% if F.missing %}x{% endif %}"):
            with self.subTest(text=text), self.assertRaises(jinja2.UndefinedError):
                sb.fill_facts(text, {"a": "1"})
        with self.assertRaisesRegex(jinja2.UndefinedError, r"^line 3: no fact 'missing'$"):
            sb.fill_facts("{{ F.a }}\n\nText {{ F.missing }}.\n", {"a": "1"})


@unittest.skipUnless((ROOT / "results" / "g2" / "g2_summary.json").exists(), "results not available")
class TestSiteBuild(unittest.TestCase):
    def test_build_refuses_a_foreign_directory(self):
        mod = builder()
        with tempfile.TemporaryDirectory() as d:
            precious = Path(d) / "keep.txt"
            precious.write_text("not a build")
            with self.assertRaises(SystemExit):
                mod.main(["--out", d])
            self.assertTrue(precious.is_file())

    def test_build_fails_on_a_missing_fact_or_reference(self):
        mod = builder()
        mod.load_facts = lambda: {"a": {"value": "1", "raw": 1, "source": "results/g2/g2_summary.json :: a", "module": "test"}}
        for body, needle in (("Value: {{ F.no_such_fact_xyz }}", "no_such_fact_xyz"), ("See [@fig-nowhere].", "fig-nowhere")):
            with self.subTest(body=body), tempfile.TemporaryDirectory() as d:
                mod.REPORT = Path(d) / "report.md"
                mod.REPORT.write_text(f"# Title\n\n{body}\n")
                with self.assertRaises(SystemExit) as cm:
                    mod.main(["--out", str(Path(d) / "site")])
                self.assertIn(needle, str(cm.exception))

    def test_build_with_facts(self):
        mod = builder()
        src = "results/g2/g2_summary.json :: arrays.opm_dense.channels"
        mod.load_facts = lambda: {"n_dense": {"value": "208 <sites>", "raw": 208, "source": src, "module": "test"},
                                  "unused": {"value": "7", "raw": 7, "source": "configs/g2.toml :: x", "module": "test"}}
        with tempfile.TemporaryDirectory() as d:
            mod.REPORT = Path(d) / "report.md"
            mod.REPORT.write_text("# Title with D~child~\n\n## 1. Results {#results}\n\nThe array has {{ F.n_dense }}.\n")
            out = Path(d) / "site"
            mod.main(["--out", str(out)])
            index = (out / "index.html").read_text()
            self.assertIn("<title>Title with Dchild · OPM vs SQUID MEG</title>", index)
            self.assertIn("<h1>Title with D<sub>child</sub></h1>", index)
            self.assertIn("<p>The array has 208 &lt;sites&gt;.</p>", index)
            numbers = (out / "numbers.html").read_text()
            self.assertIn("The report uses 1 of the 2 facts available", numbers)
            self.assertIn('<tr><td><code>n_dense</code></td><td>208 &lt;sites&gt;</td><td><a href="data/g2/g2_summary.json">'
                          'results/g2/g2_summary.json</a> :: arrays.opm_dense.channels</td></tr>', numbers)
            self.assertNotIn("unused", numbers)
            self.assertEqual(sb.check_links(out), [])

    def test_build(self):
        mod = builder()
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "site"
            mod.main(["--out", str(out)])  # raises on any broken link or anchor
            pages = {p.name for p in out.glob("*.html")}
            self.assertTrue({"index.html", "summary.html", "numbers.html", "benchmarks.html", "adult.html", "pediatric.html",
                             "epilepsy.html", "methods.html", "register.html", "reproduce.html"} <= pages)
            manifest = json.loads((out / "data" / "MANIFEST.json").read_text())
            for m in manifest:
                self.assertEqual(sb.sha256(out / m["href"]), sb.sha256(ROOT / "results" / m["path"]))
            copied = [p for p in out.rglob("*") if p.is_file() and p.name != mod.MARKER]
            # only PNG figures, result tables and the page assets: no PDF/SVG (licensed fonts), no data or caches
            self.assertEqual({p.suffix for p in copied} - {".html", ".css", ".png", ".json", ".csv", ".md"}, set())
            self.assertFalse(any(p.suffix in (".pdf", ".svg") for p in copied))
            # the former overview is the project summary
            summary = (out / "summary.html").read_text()
            self.assertIn(f"frozen as {mod.FROZEN_TAG}", summary)
            self.assertIn(mod.CURRENT_TAG, summary)
            g2 = json.loads((ROOT / "results" / "g2" / "g2_summary.json").read_text())
            self.assertIn(f"dense {g2['arrays']['opm_dense']['channels']}-site", summary)
            # the landing page is the report
            facts = mod.load_facts()
            text, used = sb.fill_facts(mod.REPORT.read_text(), {k: f["value"] for k, f in facts.items()})
            index = (out / "index.html").read_text()
            title = next(line for line in text.splitlines() if line.strip())
            self.assertIn(f"<h1>{sb.inline(title[2:].strip())}</h1>", index)
            self.assertIn('<main class="report">', index)
            for kind, start in (("Figure", "::: figure"), ("Table", "Table:")):  # numbered in order
                n = sum(line.startswith(start) for line in text.splitlines())
                self.assertEqual(re.findall(rf"<strong>{kind} (\d+)\.</strong>", index), [str(k) for k in range(1, n + 1)])
            self.assertIn('<nav class="toc" aria-label="Contents">', index)
            self.assertIn('aria-current="page">Report</a>', index)
            for ident in re.findall(r"\{#([^{}\s]+)\}", text):  # headings, figures and tables keep their ids
                self.assertIn(f'id="{ident}"', index)
            for path in re.findall(r"!\[[^\]]*\]\(results/([^)\s]+)\)", text):
                self.assertTrue((out / "figures" / path).is_file())
                self.assertIn(f'src="figures/{path}"', index)
            self.assertNotRegex(index, r"\[@|\{\{|\{%")  # every cross-reference and fact resolved
            self.assertIn("Report: Laiba Khan.", index)
            self.assertIn('Every number is read from its result files (see <a href="numbers.html">Number provenance</a>)', index)
            for page in pages:
                self.assertIn('<meta name="robots" content="noindex, nofollow">', (out / page).read_text())
            # number provenance: every fact the report uses, with its value and source
            numbers = (out / "numbers.html").read_text()
            self.assertIn(f"The report uses {len(used):,} of the {len(facts):,} facts available", numbers)
            for name in used:
                self.assertIn(f"<code>{name}</code>", numbers)
                self.assertIn(f"<td>{html.escape(str(facts[name]['value']))}</td>", numbers)
            self.assertEqual(sb.check_links(out), [])


if __name__ == "__main__":
    unittest.main()
