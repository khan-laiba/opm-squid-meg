"""The overlap audit of paper/export_figures.py, which runs on the layout of the PDF's own render as well as on the
screen's and checks texts against texts and against the symbols of the keys. Fig. 4 drawn with the row label of an
earlier revision, "as background and additional targets", is reported: its longer label column moved panel C's axis
label into a key symbol of panel D in the saved PDF, though not on the screen. The figure as drawn now is not."""
import contextlib
import importlib.util
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAME = "Figure_R17_noise_checks"
SCRIPT = ROOT / "scripts" / "report_figures_noise.py"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "paper" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@contextlib.contextmanager
def exporting(ex, out: Path):
    """The export's hooks in force (its savefig and text recording, the drawing scripts' output in out), and everything
    they change restored afterwards: the patched methods, the style module's settings, rcParams, the working directory
    and the module registry."""
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.text import Text
    style = ex.style
    saved = dict(savefig=Figure.savefig, draw=Text.draw, export=ex.EXPORT, out=style.OUT, save=style.save, dpi=style.DPI,
                 colors=dict(style.ANAT_COLOR), short=dict(style.ANAT_SHORT), cwd=os.getcwd(),
                 module=sys.modules.get("report_figures_noise"))
    ex.EXPORT = out
    os.chdir(ROOT)  # the drawing scripts read results/ relative to the repository
    try:
        with matplotlib.rc_context():
            ex.setup()
            yield
    finally:
        Figure.savefig, Text.draw, ex.EXPORT = saved["savefig"], saved["draw"], saved["export"]
        style.OUT, style.save, style.DPI = saved["out"], saved["save"], saved["dpi"]
        style.ANAT_COLOR.clear()
        style.ANAT_COLOR.update(saved["colors"])
        style.ANAT_SHORT.clear()
        style.ANAT_SHORT.update(saved["short"])
        os.chdir(saved["cwd"])
        if saved["module"] is None:
            sys.modules.pop("report_figures_noise", None)
        else:
            sys.modules["report_figures_noise"] = saved["module"]


def draw_fig4(ex, long_label: bool) -> list:
    """Fig. 4 drawn as the export draws it (draw_noise_checks()), optionally with the earlier, longer row label; the
    audit's findings for it."""
    src = SCRIPT.read_text()
    if long_label:
        assert src.count('"as both"') == 2, "the row label and its lookup"
        src = src.replace('"as both"', '"as background and additional targets"')
    mod = types.ModuleType("report_figures_noise")
    mod.__file__ = str(SCRIPT)  # its paths resolve from the repository, as when the export loads it
    sys.modules["report_figures_noise"] = mod
    exec(compile(src, mod.__file__, "exec"), mod.__dict__)
    real = ex.load
    ex.load = lambda name: mod if name == "report_figures_noise" else real(name)
    try:
        ex.draw_noise_checks()
    finally:
        ex.load = real
    return list(ex.CLASH[NAME])


class TestOverlapAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ex = _load("export_figures")

    def test_label_running_into_a_key_symbol_in_the_pdf_is_reported(self):
        with tempfile.TemporaryDirectory() as d, exporting(self.ex, Path(d)):
            found = draw_fig4(self.ex, long_label=True)
        self.assertTrue(any("a key's symbol" in c for c in found), found)

    def test_figure_as_drawn_now_has_no_overlap(self):
        with tempfile.TemporaryDirectory() as d, exporting(self.ex, Path(d)):
            self.assertEqual(draw_fig4(self.ex, long_label=False), [])


if __name__ == "__main__":
    unittest.main()
