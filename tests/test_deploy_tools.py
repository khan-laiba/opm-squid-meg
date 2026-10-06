"""Deployment tools, offline: the gh-pages commit of scripts/deploy_pages.py on a temporary git repository and
scripts/check_live_site.py on a site served by a local HTTP server."""
import contextlib
import functools
import hashlib
import importlib.util
import io
import os
import subprocess
import tempfile
import threading
import unittest
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


deploy = _load("deploy_pages")
live = _load("check_live_site")


class TestDeployPages(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        (self.tmp / "gitconfig").write_text("[user]\n\tname = Global Person\n\temail = global@example.invalid\n")
        # a global config of its own: none of the user's signing, hooks or excludes, and an identity never to be used
        env = mock.patch.dict(os.environ, {"GIT_CONFIG_GLOBAL": str(self.tmp / "gitconfig"), "GIT_CONFIG_NOSYSTEM": "1"})
        env.start()
        self.addCleanup(env.stop)
        self.repo = self.tmp / "repo"
        subprocess.run(["git", "init", "-q", "-b", "main", str(self.repo)], check=True)
        self.git("config", "user.name", "Test Author")
        self.git("config", "user.email", "author@example.invalid")
        self.write(self.repo, {"README.md": "readme\n", "scripts/a.py": "print(1)\n"})
        self.git("add", ".")
        self.git("commit", "-q", "-m", "init")
        # as in the project's .git/info/exclude: ignore rules must not drop files of the build
        self.write(self.repo, {".git/info/exclude": "/data\n*.png\n"})
        self.build = self.tmp / "build"
        self.write(self.build, {"index.html": "<h1>v1</h1>", "data/table.json": "{}", "figures/fig.png": "png"})

    def git(self, *args):
        return subprocess.run(["git", "-C", str(self.repo), *args], check=True, capture_output=True, text=True).stdout.strip()

    @staticmethod
    def write(base, files):
        for rel, text in files.items():
            (base / rel).parent.mkdir(parents=True, exist_ok=True)
            (base / rel).write_text(text)

    def test_commit_creates_the_branch_then_a_child_without_touching_the_checkout(self):
        head, index = self.git("rev-parse", "HEAD"), hashlib.sha256((self.repo / ".git" / "index").read_bytes()).digest()
        c1 = deploy.commit_pages(self.repo, self.build, "Pages: first")
        self.assertEqual(self.git("rev-parse", "refs/heads/gh-pages"), c1)
        self.assertEqual(sorted(deploy.tree_files(self.repo, c1)), [".nojekyll", "data/table.json", "figures/fig.png", "index.html"])
        self.assertEqual(self.git("cat-file", "-s", "gh-pages:.nojekyll"), "0")
        self.assertEqual(self.git("show", "gh-pages:index.html"), "<h1>v1</h1>")
        self.assertEqual(self.git("rev-list", "--count", "gh-pages"), "1")  # a root commit
        ident = "Test Author <author@example.invalid>"
        self.assertEqual(self.git("log", "-1", "--format=%an <%ae>|%cn <%ce>|%s", "gh-pages"), f"{ident}|{ident}|Pages: first")
        # the next build: one file changed, one removed, one added; the identity comes from the git config
        (self.build / "figures" / "fig.png").unlink()
        self.write(self.build, {"index.html": "<h1>v2</h1>", "new.html": "<p>new</p>"})
        with mock.patch.dict(os.environ, {"GIT_AUTHOR_NAME": "Someone Else", "GIT_COMMITTER_EMAIL": "else@example.invalid"}):
            c2 = deploy.commit_pages(self.repo, self.build, "Pages: second")
        self.assertEqual(self.git("rev-parse", "gh-pages"), c2)
        self.assertEqual(self.git("rev-parse", "gh-pages^"), c1)
        self.assertEqual(sorted(deploy.tree_files(self.repo, c2)), [".nojekyll", "data/table.json", "index.html", "new.html"])
        self.assertEqual(self.git("show", "gh-pages:index.html"), "<h1>v2</h1>")
        self.assertEqual(self.git("log", "-1", "--format=%an <%ae>|%cn <%ce>", "gh-pages"), f"{ident}|{ident}")
        # the checkout is untouched: the same branch, commit and index file, and a clean working tree
        self.assertEqual(hashlib.sha256((self.repo / ".git" / "index").read_bytes()).digest(), index)
        self.assertEqual(self.git("symbolic-ref", "HEAD"), "refs/heads/main")
        self.assertEqual(self.git("rev-parse", "HEAD"), head)
        self.assertEqual(self.git("status", "--porcelain", "--untracked-files=all"), "")

    def test_parent_from_origin_without_a_local_branch(self):
        c1 = deploy.commit_pages(self.repo, self.build, "Pages: first")
        self.git("update-ref", "refs/remotes/origin/gh-pages", c1)
        self.git("update-ref", "-d", "refs/heads/gh-pages")
        deploy.commit_pages(self.repo, self.build, "Pages: second")
        self.assertEqual(self.git("rev-parse", "gh-pages^"), c1)

    def test_identity_only_from_the_repository_config(self):
        self.git("config", "--unset", "user.email")
        with self.assertRaises(RuntimeError):  # the global identity is not used
            deploy.commit_pages(self.repo, self.build, "Pages: first")
        self.assertEqual(deploy.rev(self.repo, "refs/heads/gh-pages"), "")

    def test_refuses_a_checked_out_branch(self):
        c1 = deploy.commit_pages(self.repo, self.build, "Pages: first")
        self.git("checkout", "-q", "gh-pages")
        self.write(self.build, {"index.html": "<h1>v2</h1>"})
        with self.assertRaises(RuntimeError):
            deploy.commit_pages(self.repo, self.build, "Pages: second")
        self.assertEqual(self.git("rev-parse", "gh-pages"), c1)

    def test_source_commit_marks_uncommitted_sources(self):
        short = self.git("rev-parse", "--short", "HEAD")
        self.assertEqual(deploy.source_commit(self.repo), (short, []))
        self.write(self.repo, {"notes/x.txt": "a = 1\n"})  # outside the source paths
        self.assertEqual(deploy.source_commit(self.repo), (short, []))
        self.write(self.repo, {"docs/new.md": "new\n"})  # untracked, in a source path
        self.assertEqual(deploy.source_commit(self.repo), (short + "+dirty", ["?? docs/new.md"]))
        (self.repo / "docs" / "new.md").unlink()
        self.write(self.repo, {"scripts/a.py": "print(2)\n"})
        self.assertEqual(deploy.source_commit(self.repo), (short + "+dirty", [" M scripts/a.py"]))


class _Handler(SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/forbidden":
            self.send_error(403)
        else:
            super().do_GET()

    def log_message(self, *args):
        pass


class TestCheckLiveSite(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(_Handler, directory=cls.tmp.name))
        threading.Thread(target=cls.server.serve_forever, kwargs=dict(poll_interval=0.05), daemon=True).start()
        port = cls.server.server_address[1]
        cls.base, cls.ext = f"http://127.0.0.1:{port}", f"http://localhost:{port}"  # same server, another host: external
        site = {
            "outside.html": "<p>outside the site's path prefix</p>",
            "site/index.html": ('<link rel="stylesheet" href="static/style.css"><h1 id="local">Report</h1>'
                                '<a href="b.html">b</a> <a href="b.html#sec">good anchor</a> <a href="b.html#nope">missing anchor</a> '
                                '<a href="#local">this page</a> <img src="img/ok.png"> <img src="img/missing.png"> '
                                '<a href="data/table.csv">table</a> <a href="../outside.html">outside</a> '
                                '<a href="https://example.org/paper">paper</a> <a href="mailto:someone@example.org">mail</a>'),
            "site/b.html": '<h2 id="sec">Section</h2><a href="index.html#local">back</a>',
            "site/static/style.css": "body {}",
            "site/img/ok.png": "png",
            "site/data/table.csv": "a,b\n1,2\n",
            "ext/index.html": (f'<a href="{cls.ext}/site/img/ok.png">ok</a> <a href="{cls.ext}/site/img">redirect</a> '
                               f'<a href="{cls.ext}/forbidden">blocked</a> <a href="{cls.ext}/gone">gone</a>'),
            "blocked/index.html": f'<a href="{cls.ext}/forbidden">blocked only</a>',
        }
        for rel, text in site.items():
            (Path(cls.tmp.name) / rel).parent.mkdir(parents=True, exist_ok=True)
            (Path(cls.tmp.name) / rel).write_text(text)
        cls.patches = [mock.patch.object(live, "RETRY_DELAY", 0), mock.patch.dict(os.environ, {"no_proxy": "*", "NO_PROXY": "*"})]
        for p in cls.patches:
            p.start()

    @classmethod
    def tearDownClass(cls):
        for p in cls.patches:
            p.stop()
        cls.server.shutdown()
        cls.server.server_close()
        cls.tmp.cleanup()

    def main(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            status = live.main(list(argv))
        return status, out.getvalue()

    def test_internal_failures_found_and_good_links_pass(self):
        rep = live.check(f"{self.base}/site/")
        failures = dict(sorted(rep["failures"]))
        self.assertEqual(list(failures), [f"404 {self.base}/site/img/missing.png (File not found)",
                                          f"missing anchor {self.base}/site/b.html#nope"])
        # both linked from the start page, reached as site/ and (from b.html) as site/index.html
        self.assertEqual(set().union(*failures.values()), {f"{self.base}/site/", f"{self.base}/site/index.html"})
        # pages: site/ (= site/index.html) and b.html; assets: the style sheet, two images, the table;
        # anchors: b.html#sec, b.html#nope, site/#local and site/index.html#local
        self.assertEqual((rep["pages"], rep["assets"], rep["anchors"]), (2, 4, 4))
        self.assertEqual(rep["external"].keys(), {f"{self.base}/outside.html", "https://example.org/paper"})
        self.assertTrue(all(v["verdict"] is None for v in rep["external"].values()))  # not fetched without --external
        self.assertEqual(live.check(f"{self.base}/site"), rep)  # the same after the server's redirect to site/

    def test_exit_status_and_summary(self):
        status, out = self.main(f"{self.base}/site/")
        self.assertEqual(status, 1)
        self.assertIn("2 pages, 4 internal assets, 4 anchors checked, 2 internal failures", out)
        self.assertIn("external links: 2, not checked", out)

    def test_external_links(self):
        with mock.patch.object(live, "PUBLISHERS", live.PUBLISHERS + ("localhost",)):
            rep = live.check(f"{self.base}/ext/", external=True)
            self.assertEqual(self.main(f"{self.base}/blocked/", "--external")[0], 0)  # blocked links do not fail
            status, out = self.main(f"{self.base}/ext/", "--external")
        self.assertEqual({u: v["verdict"] for u, v in rep["external"].items()},
                         {f"{self.ext}/site/img/ok.png": "ok", f"{self.ext}/site/img": "ok",  # the second after a redirect
                          f"{self.ext}/forbidden": "blocked", f"{self.ext}/gone": "failed"})
        self.assertEqual(rep["failures"], [])
        self.assertEqual(status, 1)
        self.assertIn("external links: 4; 2 ok, 1 blocked to scripts (check in a browser), 1 failed", out)
        # a 403 from a host that is not a publisher fails
        self.assertEqual(self.main(f"{self.base}/blocked/", "--external")[0], 1)

    def test_classify(self):
        c = live.classify
        self.assertEqual(c(200, "https://example.org/a", "https://example.org/a"), "ok")
        self.assertEqual(c(302, "https://example.org/a", "https://example.org/b"), "ok")
        self.assertEqual(c(403, "https://doi.org/10.1000/x", "https://onlinelibrary.wiley.com/doi/x"), "blocked")
        self.assertEqual(c(429, "https://www.biorxiv.org/content/x", "https://www.biorxiv.org/content/x"), "blocked")
        self.assertEqual(c(401, "https://example.org/a", "https://www.sciencedirect.com/x"), "blocked")  # after a redirect
        self.assertEqual(c(403, "https://example.org/a", "https://example.org/a"), "failed")
        self.assertEqual(c(403, "https://notdoi.org/x", "https://notdoi.org/x"), "failed")
        self.assertEqual(c(404, "https://doi.org/10.1000/x", "https://doi.org/10.1000/x"), "failed")
        self.assertEqual(c(500, "https://doi.org/10.1000/x", "https://doi.org/10.1000/x"), "failed")
        self.assertEqual(c(None, "https://doi.org/10.1000/x", "https://doi.org/10.1000/x"), "failed")


if __name__ == "__main__":
    unittest.main()
