import contextlib
import html
import importlib.util
import io
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

    def test_comments_are_dropped_and_text_after_them_is_an_error(self):
        h = sb.md_to_html("Before.\n<!-- one line -->\nAfter one.\n\n<!-- two\nlines -->\n\nAfter two.\n  <!-- indented -->  \nEnd.\n")
        self.assertEqual(h, "<p>Before.</p>\n<p>After one.</p>\n<p>After two.</p>\n<p>End.</p>")
        # text after '-->' on the comment's last line used to vanish with the comment: now an error naming the line
        for md, line in (("Text.\n\n<!-- B3: a note --> The sentence that would vanish.\n", 3),
                         ("# T\n\n<!-- a note\n   over two lines --> lost\n\nNext.\n", 4)):
            with self.subTest(md=md), self.assertRaisesRegex(ValueError, f"^line {line}: text after '-->' would be dropped"):
                sb.md_to_html(md)
        with self.assertRaisesRegex(ValueError, "^line 3: comment not closed by '-->'"):  # the rest of the document is not lost
            sb.md_to_html("Text.\n\n<!-- never closed\n\nA paragraph.\n")

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

    def test_numbers_with_a_prefix_and_references_between_documents(self):
        report = ("See [@fig-s-qc], [@tab-s-boot] and [@fig-a].\n\n"
                  "::: figure {#fig-a}\n![x](results/g9/a.png)\n:::\n")
        supplement = ("As in [@fig-a]; see [@fig-s-qc] and `[@fig-a]`.\n\n"
                      "::: figure {#fig-s-qc}\n![y](results/g9/b.png)\nCaption, see [@tab-s-boot].\n:::\n\n"
                      "Table: Resamples {#tab-s-boot}\n| a |\n|---|\n| 1 |\n\n"
                      "::: figure {#fig-s-two}\n![z](results/g9/b.png)\n:::\n")
        rh, rlab = sb.md_blocks(report, image=lambda p: sb.publish_png(p, self.res, self.out))
        shtml, slab = sb.md_blocks(supplement, image=lambda p: sb.publish_png(p, self.res, self.out), prefix="S")
        self.assertEqual(rlab, {"fig-a": "Figure 1"})
        self.assertEqual(list(slab.items()), [("fig-s-qc", "Figure S1"), ("tab-s-boot", "Table S1"), ("fig-s-two", "Figure S2")])
        self.assertIn("[@fig-s-qc]", rh)  # left for resolve_xrefs
        self.assertIn("<figcaption><strong>Figure S1.</strong> Caption, see [@tab-s-boot].</figcaption>", shtml)
        self.assertIn("<caption><strong>Table S1.</strong> Resamples</caption>", shtml)
        self.assertIn("<figcaption><strong>Figure S2.</strong></figcaption>", shtml)
        h = sb.resolve_xrefs(rh, rlab, {i: ("supplement.html", lab) for i, lab in slab.items()})
        self.assertIn('See <a href="supplement.html#fig-s-qc">Figure S1</a>, <a href="supplement.html#tab-s-boot">Table S1</a> '
                      'and <a href="#fig-a">Figure 1</a>.', h)
        s = sb.resolve_xrefs(shtml, slab, {i: ("index.html", lab) for i, lab in rlab.items()})
        self.assertIn('As in <a href="index.html#fig-a">Figure 1</a>; see <a href="#fig-s-qc">Figure S1</a> and '
                      '<code>[@fig-a]</code>.', s)  # never inside code
        self.assertIn('Caption, see <a href="#tab-s-boot">Table S1</a>.', s)
        # one document: md_to_html with a prefix and with another document's labels
        self.assertIn("<strong>Table S1.</strong>", sb.md_to_html("Table: T {#tab-x}\n| a |\n|---|\n| 1 |\n", prefix="S"))
        self.assertIn('<a href="index.html#fig-a">Figure 1</a>',
                      sb.md_to_html("See [@fig-a].", refs={"fig-a": ("index.html", "Figure 1")}))
        with self.assertRaisesRegex(ValueError, "unknown figure or table: \\[@fig-s-missing\\]"):
            sb.resolve_xrefs("See [@fig-s-missing].", rlab, {i: ("supplement.html", lab) for i, lab in slab.items()})
        with self.assertRaisesRegex(ValueError, "defined in both documents: fig-a"):  # an id may not be defined in both
            sb.resolve_xrefs(rh, rlab, {"fig-a": ("supplement.html", "Figure S3")})

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


# milestone codes and internal keys that no page title may show
INTERNAL = re.compile(r"\b(G[0-5][A-C]?|REPRO|ADAPT|NEW)\b|squid/|opm_dense|opm_matched|intrinsic")
FIG_R1 = "::: figure {#fig-r}\n![Depth](results/g2/Figure_G2_depth.png)\nCaption, see [@fig-s-qc].\n:::\n"
FIG_S1 = "::: figure {#fig-s-qc}\n![QC](results/g2/Figure_G2_bands.png)\nCaption, see [@fig-r].\n:::\n"
TAB_S1 = "Table: Resamples {#tab-s-boot}\n\n| a | b |\n|---|---|\n| 1 | 2 |\n"


def nav_titles(page: str) -> list:
    """The supplementary pages in a page's navigation, in order."""
    nav = page[page.index('<nav aria-label="Sections">'):page.index("</nav>")]
    return re.findall(r'<a [^>]*href="([^"#]+)"[^>]*>([^<]+)</a>', nav)


@unittest.skipUnless((ROOT / "results" / "g2" / "g2_summary.json").exists(), "results not available")
class TestSiteBuild(unittest.TestCase):
    @staticmethod
    def documents(mod, d, report, supplement=None):
        """A stub report and, if given, a stub supplementary text in the test's own folder (never in report/)."""
        mod.REPORT, mod.SUPPLEMENT_MD = Path(d) / "report.md", Path(d) / "supplement.md"
        mod.REPORT.write_text(report)
        if supplement is not None:
            mod.SUPPLEMENT_MD.write_text(supplement)

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
        s1 = "# S1. Supplementary text\n\n"
        cases = (("Value: {{ F.no_such_fact_xyz }}", None, "report.md", "no_such_fact_xyz"),
                 ("See [@fig-nowhere].", None, "report.md", "fig-nowhere"),
                 ("See [@fig-s-qc].", None, "report.md", "fig-s-qc"),  # a supplementary figure, but no supplementary text
                 ("Text.", s1 + "Value: {{ F.no_such_fact_s }}\n", "supplement.md", "no_such_fact_s"),
                 ("Text.", s1 + "See [@fig-nowhere].\n", "supplement.md", "fig-nowhere"),
                 ("Text.", "No title line.\n", "supplement.md", "title"),
                 # a comment that would swallow text: the error names the file and the line in it
                 ("<!-- B3: a note --> The sentence that would vanish.", None, "report.md", "line 3: text after '-->'"),
                 ("Text.", s1 + "<!-- a\nnote --> lost\n", "supplement.md", "line 4: text after '-->'"),
                 (FIG_R1.replace("fig-s-qc", "fig-r"), s1 + FIG_R1.replace("fig-s-qc", "fig-r"), "report.md",
                  "defined in both documents: fig-r"))
        for body, supplement, where, needle in cases:
            with self.subTest(body=body, supplement=supplement), tempfile.TemporaryDirectory() as d:
                self.documents(mod, d, f"# Title\n\n{body}\n", supplement)
                with self.assertRaises(SystemExit) as cm, contextlib.redirect_stdout(io.StringIO()):
                    mod.main(["--out", str(Path(d) / "site")])
                self.assertIn(needle, str(cm.exception))
                self.assertIn(where, str(cm.exception))

    def test_build_with_facts(self):
        mod = builder()
        src = "results/g2/g2_summary.json :: arrays.opm_dense.channels"
        mod.load_facts = lambda: {"n_dense": {"value": "208 <sites>", "raw": 208, "source": src, "module": "test"},
                                  "unused": {"value": "7", "raw": 7, "source": "configs/g2.toml :: x", "module": "test"}}
        with tempfile.TemporaryDirectory() as d:
            self.documents(mod, d, "# Title with D~child~\n\n## 1. Results {#results}\n\nThe array has {{ F.n_dense }}.\n")
            out = Path(d) / "site"
            with contextlib.redirect_stdout(io.StringIO()) as log:
                mod.main(["--out", str(out)])
            # without report/supplement.md: a warning, no S1 page and no S1 in the navigation
            self.assertIn(f"warning: {mod.SUPPLEMENT_MD} not found; building without the supplementary text (S1)", log.getvalue())
            self.assertFalse((out / "supplement.html").exists())
            index = (out / "index.html").read_text()
            self.assertEqual(nav_titles(index)[1][1], "S2 Adult reference benchmarks")
            self.assertNotIn("supplement.html", index)
            self.assertIn("<title>Title with Dchild · OPM vs SQUID MEG</title>", index)
            self.assertIn("<h1>Title with D<sub>child</sub></h1>", index)
            self.assertIn("<p>The array has 208 &lt;sites&gt;.</p>", index)
            numbers = (out / "numbers.html").read_text()
            self.assertIn("The report uses 1 of the 2 facts available", numbers)
            self.assertIn('<tr><td><code>n_dense</code></td><td>208 &lt;sites&gt;</td><td><a href="data/g2/g2_summary.json">'
                          'results/g2/g2_summary.json</a> :: arrays.opm_dense.channels</td></tr>', numbers)
            self.assertNotIn("unused", numbers)
            self.assertEqual(sb.check_links(out), [])

    def test_build_with_supplementary_text(self):
        mod = builder()
        g2 = "results/g2/g2_summary.json :: arrays."
        mod.load_facts = lambda: {
            "n_dense": {"value": "208", "raw": 208, "source": g2 + "opm_dense.channels", "module": "test"},
            "n_matched": {"value": "98", "raw": 98, "source": g2 + "opm_matched.channels", "module": "test"},
            "unused": {"value": "7", "raw": 7, "source": "configs/g2.toml :: x", "module": "test"}}
        report = ("# Title\n\n## 1. Results {#results}\n\nThe dense array has {{ F.n_dense }} sites ([@fig-r]; [@fig-s-qc], "
                  "[@tab-s-boot]).\n\n" + FIG_R1)
        supplement = ("# S1. Supplementary text\n\n## A. Conventions {#s-conventions}\n\nMatched {{ F.n_matched }}, dense "
                      "{{ F.n_dense }} ([@fig-r], [@fig-s-qc]).\n\n" + FIG_S1 + "\n" + TAB_S1 + "\n## B. More\n\nText.\n")
        with tempfile.TemporaryDirectory() as d:
            self.documents(mod, d, report, supplement)
            out = Path(d) / "site"
            with contextlib.redirect_stdout(io.StringIO()) as log:
                mod.main(["--out", str(out)])
            self.assertNotIn("warning", log.getvalue())
            index, supp = (out / "index.html").read_text(), (out / "supplement.html").read_text()
            # the supplementary text: its own title, the report's layout and contents, figures and tables numbered S1, S2, ...
            self.assertIn("<title>S1. Supplementary text · OPM vs SQUID MEG</title>", supp)
            self.assertIn('<main class="report">\n<h1>S1. Supplementary text</h1>', supp)
            self.assertIn('<nav class="toc" aria-label="Contents">', supp)
            self.assertIn('<h2 id="s-conventions">A. Conventions</h2>', supp)
            self.assertIn("<p>Matched 98, dense 208 (", supp)
            self.assertEqual(re.findall(r"<strong>(Figure|Table) (S?\d+)\.</strong>", supp), [("Figure", "S1"), ("Table", "S1")])
            self.assertEqual(re.findall(r"<strong>(Figure|Table) (S?\d+)\.</strong>", index), [("Figure", "1")])
            self.assertTrue((out / "figures" / "g2" / "Figure_G2_bands.png").is_file())
            # cross-references within and between the two documents
            self.assertIn('(<a href="#fig-r">Figure 1</a>; <a href="supplement.html#fig-s-qc">Figure S1</a>, '
                          '<a href="supplement.html#tab-s-boot">Table S1</a>)', index)
            self.assertIn('Caption, see <a href="supplement.html#fig-s-qc">Figure S1</a>.', index)
            self.assertIn('(<a href="index.html#fig-r">Figure 1</a>, <a href="#fig-s-qc">Figure S1</a>)', supp)
            self.assertIn('Caption, see <a href="index.html#fig-r">Figure 1</a>.', supp)
            # navigation and page titles: S1-S9, then the project summary
            expected = [("supplement.html", "S1 Supplementary text"), ("benchmarks.html", "S2 Adult reference benchmarks"),
                        ("adult.html", "S3 Realistic adult comparison"), ("pediatric.html", "S4 Smaller heads"),
                        ("epilepsy.html", "S5 Simulated interictal spikes and head motion"),
                        ("methods.html", "S6 Methods, uncertainty and limitations"),
                        ("register.html", "S7 Parameters, provenance and assumptions"), ("numbers.html", "S8 Number provenance"),
                        ("reproduce.html", "S9 Reproduce and download"), ("summary.html", "Project summary by milestone")]
            self.assertEqual(nav_titles(index), [("index.html", "Report")] + expected)
            self.assertIn('<a href="supplement.html" aria-current="page">S1 Supplementary text</a>', supp)
            for page, label in expected[1:]:
                number, _, name = label.partition(" ")
                title = f"{number}. {name}" if number.startswith("S") else label
                self.assertIn(f"<h1>{html.escape(title)}</h1>", (out / page).read_text())
            # number provenance: the facts of both documents, with the documents that use each
            numbers = (out / "numbers.html").read_text()
            self.assertIn("The report uses 1 and the supplementary text 2 of the 3 facts available (2 different facts, 1 of them in "
                          "both)", numbers)
            self.assertIn('<a href="supplement.html">supplementary text</a> (S1)', numbers)
            self.assertIn('<tr><td><code>n_dense</code></td><td>208</td><td><a href="data/g2/g2_summary.json">results/g2/'
                          'g2_summary.json</a> :: arrays.opm_dense.channels</td><td>report, S1</td></tr>', numbers)
            self.assertIn("</a> :: arrays.opm_matched.channels</td><td>S1</td></tr>", numbers)
            self.assertLess(numbers.index("<code>n_dense</code>"), numbers.index("<code>n_matched</code>"))  # the report's first
            self.assertNotIn("<code>unused</code>", numbers)
            self.assertEqual(sb.check_links(out), [])

    def test_build(self):
        mod = builder()
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "site"
            with contextlib.redirect_stdout(io.StringIO()):
                mod.main(["--out", str(out)])  # raises on any broken link or anchor
            pages = {p.name for p in out.glob("*.html")}
            has_supplement = mod.SUPPLEMENT_MD.is_file()
            self.assertEqual("supplement.html" in pages, has_supplement)
            self.assertTrue({"index.html", "summary.html", "numbers.html", "benchmarks.html", "adult.html", "pediatric.html",
                             "epilepsy.html", "methods.html", "register.html", "reproduce.html"} <= pages)
            manifest = json.loads((out / "data" / "MANIFEST.json").read_text())
            for m in manifest:
                self.assertEqual(sb.sha256(out / m["href"]), sb.sha256(ROOT / "results" / m["path"]))
            # a table records the commit of the run that wrote it (its first line), as on the download page
            for m in manifest:
                if m["path"].endswith(".csv"):
                    first = (ROOT / "results" / m["path"]).read_text().split("\n", 1)[0]
                    found = re.search(r"\|\s*commit\s+(\S+)\s*$", first)
                    if found:
                        self.assertEqual(m["commit"], found.group(1), m["path"])
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
            values = {k: f["value"] for k, f in facts.items()}
            text, used = sb.fill_facts(mod.REPORT.read_text(), values)
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
            # the supplementary text (S1), if there is one: numbered S1, S2, ..., every reference and fact resolved
            used_s = []
            if has_supplement:
                text_s, used_s = sb.fill_facts(mod.SUPPLEMENT_MD.read_text(), values)
                supp = (out / "supplement.html").read_text()
                self.assertIn('<main class="report">', supp)
                for kind, start in (("Figure", "::: figure"), ("Table", "Table:")):
                    n = sum(line.startswith(start) for line in text_s.splitlines())
                    self.assertEqual(re.findall(rf"<strong>{kind} (S?\d+)\.</strong>", supp), [f"S{k}" for k in range(1, n + 1)])
                for ident in re.findall(r"\{#([^{}\s]+)\}", text_s):
                    self.assertIn(f'id="{ident}"', supp)
                self.assertNotRegex(supp, r"\[@|\{\{|\{%")
            # navigation: the supplementary pages in order, S1 (if built) to S9, then the project summary
            built = [(f, f"{n} {t}" if n else t) for f, n, t in mod.SUPPLEMENT if f in pages]
            self.assertEqual(nav_titles(index), [("index.html", "Report")] + built)
            for page in pages:
                text_p = (out / page).read_text()
                self.assertIn('<meta name="robots" content="noindex, nofollow">', text_p)
                for t in re.findall(r"<title>(.*?)</title>|<h1>(.*?)</h1>", text_p):  # page titles in the manuscript's terms
                    self.assertNotRegex("".join(t), INTERNAL, page)
            for f, n, t in mod.SUPPLEMENT:
                if f in pages and f != "supplement.html":
                    self.assertIn(f"<h1>{html.escape(mod.page_title(f))}</h1>", (out / f).read_text())
            # number provenance: every fact the report and S1 use, with its value and source
            numbers = (out / "numbers.html").read_text()
            if has_supplement:
                self.assertIn(f"The report uses {len(used):,} and the supplementary text {len(used_s):,} of the "
                              f"{len(facts):,} facts available", numbers)
            else:
                self.assertIn(f"The report uses {len(used):,} of the {len(facts):,} facts available", numbers)
            for name in dict.fromkeys(used + used_s):
                self.assertIn(f"<code>{name}</code>", numbers)
                self.assertIn(f"<td>{html.escape(str(facts[name]['value']))}</td>", numbers)
            self.assertEqual(sb.check_links(out), [])


if __name__ == "__main__":
    unittest.main()
