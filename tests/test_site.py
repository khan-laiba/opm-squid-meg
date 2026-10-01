import json
import tempfile
import unittest
from pathlib import Path

from opmsquid import sitebuild as sb

ROOT = Path(__file__).resolve().parents[1]


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


@unittest.skipUnless((ROOT / "results" / "g2" / "g2_summary.json").exists(), "results not available")
class TestSiteBuild(unittest.TestCase):
    def test_build_refuses_a_foreign_directory(self):
        import importlib.util

        spec = importlib.util.spec_from_file_location("build_site", ROOT / "scripts" / "build_site.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with tempfile.TemporaryDirectory() as d:
            precious = Path(d) / "keep.txt"
            precious.write_text("not a build")
            with self.assertRaises(SystemExit):
                mod.main(["--out", d])
            self.assertTrue(precious.is_file())

    def test_build(self):
        import importlib.util

        spec = importlib.util.spec_from_file_location("build_site", ROOT / "scripts" / "build_site.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "site"
            mod.main(["--out", str(out)])  # raises on any broken link or anchor
            pages = {p.name for p in out.glob("*.html")}
            self.assertTrue({"index.html", "benchmarks.html", "adult.html", "pediatric.html", "epilepsy.html", "methods.html",
                             "register.html", "reproduce.html"} <= pages)
            manifest = json.loads((out / "data" / "MANIFEST.json").read_text())
            for m in manifest:
                self.assertEqual(sb.sha256(out / m["href"]), sb.sha256(ROOT / "results" / m["path"]))
            copied = [p for p in out.rglob("*") if p.is_file() and p.name != mod.MARKER]
            # only PNG figures, result tables and the page assets: no PDF/SVG (licensed fonts), no data or caches
            self.assertEqual({p.suffix for p in copied} - {".html", ".css", ".png", ".json", ".csv", ".md"}, set())
            self.assertFalse(any(p.suffix in (".pdf", ".svg") for p in copied))
            self.assertIn(f"frozen as {mod.FROZEN_TAG}", (out / "index.html").read_text())
            g2 = json.loads((ROOT / "results" / "g2" / "g2_summary.json").read_text())
            self.assertIn(f"dense {g2['arrays']['opm_dense']['channels']}-site", (out / "index.html").read_text())


if __name__ == "__main__":
    unittest.main()
