"""The HTML converter paper/build_html.py on small rendered sources: a table continued with \\ContinuedFloat keeps
its number, notes set by hand as \\tablenotes sets them stay with their table, and a line break in a cell stays one
in the page's text."""
import importlib.util
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "paper" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # its dataclasses look their module up there
    spec.loader.exec_module(mod)
    return mod


H = _load("build_html")


def table(env, head, notes=""):
    return (f"\\begin{{{env}}}\n{head}\n\\begin{{tabular}}{{L{{2cm}}L{{2cm}}}}\n\\toprule\nHead & Value \\\\\n"
            f"\\midrule\nAdult & +1.00\\newline [+0.82, +1.17] \\\\\n\\bottomrule\n\\end{{tabular}}\n{notes}\n"
            f"\\end{{{env}}}\n\n")


def render(body, labels):
    tex = "\\renewcommand{\\thetable}{S\\arabic{table}}\n\\begin{document}\n" + body + "\n\\end{document}\n"
    r = H.Renderer("supplementary", tex, (labels, {}), [], Path("."), Path("."))
    return r, r.render()


def figures(page):
    """(id, label, caption, [notes]) of every table of a page."""
    return [(m.group(1), m.group(2), m.group(3), re.findall(r'<div class="(tablenotes[^"]*)">(.*?)</div>', m.group(4)))
            for m in re.finditer(r'<figure class="table[^"]*" id="([^"]+)"><figcaption><span class="label">([^<]+)'
                                 r'</span> ([^<]*)</figcaption>(.*?)</figure>', page, re.S)]


class TestTables(unittest.TestCase):
    def test_continued_table_keeps_its_number_and_its_notes(self):
        body = (table("sidewaystable", "\\caption{Studies.}\\label{tab:S1}\n\\scriptsize",
                      "\\par\\smallskip\n{\\footnotesize\\raggedleft Continued on the next page.\\par}")
                + table("sidewaystable", "\\ContinuedFloat\n\\caption[]{(Continued.)}",
                        "\\par\\smallskip\n{\\footnotesize\\raggedright Direction: as reported. "
                        "$^{\\mathrm{d}}$Derived; SD, standard deviation.\\par}")
                + table("table", "\\caption{Seeds.}\\label{tab:S2}", "\\tablenotes{Seeds as drawn.}")
                + "Table~\\ref{tab:S2} follows Table~\\ref{tab:S1}.\n")
        r, page = render(body, {"tab:S1": "S1", "tab:S2": "S2"})
        self.assertEqual(figures(page), [
            ("tab-S1", "Table\u00a0S1.", "Studies.", []),
            ("tab-S1-continued", "Table\u00a0S1.", "(Continued.)",
             [("tablenotes", "Direction: as reported. <sup>d</sup>Derived; SD, standard deviation.")]),
            ("tab-S2", "Table\u00a0S2.", "Seeds.", [("tablenotes", "Seeds as drawn.")])])
        self.assertIn('<a class="xref" href="#tab-S2">S2</a> follows', page)
        self.assertEqual(r.warnings, [])  # every number as the .aux has it

    def test_notes_after_a_longtable(self):
        body = ("{\\scriptsize\n\\begin{longtable}{L{2cm}L{2cm}}\n\\caption{Estimands.}\\label{tab:S1}\\\\\n"
                "\\toprule\nA & B \\\\\n\\midrule\n\\endfirsthead\n\\bottomrule\n\\endlastfoot\nx & y \\\\\n"
                "\\end{longtable}\n}\n\\vspace{-\\baselineskip}\n\\tablenotes{First note.}\n"
                "{\\footnotesize\\raggedright Second note.\\par}\n\nText after the table.\n")
        r, page = render(body, {"tab:S1": "S1"})
        self.assertEqual(figures(page), [("tab-S1", "Table\u00a0S1.", "Estimands.",
                                          [("tablenotes", "First note."), ("tablenotes", "Second note.")])])
        self.assertIn("<p>Text after the table.</p>", page)
        self.assertEqual(r.warnings, [])

    def test_line_break_in_a_cell(self):
        _, page = render(table("table", "\\caption{Heads.}\\label{tab:S1}"), {"tab:S1": "S1"})
        self.assertIn("<td>+1.00<br>\n[+0.82, +1.17]</td>", page)
        self.assertIn("+1.00\n[+0.82, +1.17]", H.strip_tags(page))  # the text of the page keeps the break


if __name__ == "__main__":
    unittest.main()
