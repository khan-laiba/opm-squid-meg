#!/usr/bin/env python3
"""Deploy the report to GitHub Pages: commit a fresh build as the next commit of the branch gh-pages.

scripts/build_site.py builds the site into a new temporary directory; the deployment stops if the build
fails or has a broken local link or anchor (checked again here). The build, without the builder's marker
file and with an empty .nojekyll, is committed through a temporary index (git add --work-tree, write-tree,
commit-tree, update-ref), so the working tree, the index and the current branch are never touched. The
parent is the current gh-pages commit (else origin/gh-pages, else none); author and committer are user.name
and user.email of the repository's own config (.git/config; a global identity is never used); the message
names the source commit (`git rev-parse --short HEAD`, with '+dirty' if SOURCE_PATHS have uncommitted or
untracked changes, refused unless --allow-dirty).
Nothing leaves the machine unless --push is given: then `git push origin gh-pages` (a normal push, never
forced).

Usage: .venv/bin/python scripts/deploy_pages.py [--allow-dirty] [--push]
Check: .venv/bin/python scripts/check_live_site.py https://khan-laiba.github.io/opm-squid-meg/ [--external]
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from opmsquid import sitebuild as sb  # noqa: E402

BRANCH = "gh-pages"
MARKER = ".opmsquid_site_build"  # build_site.py's mark of its own output directory (not published)
# what the build reads (README.md: its Reproduce section); a change here makes the source commit '+dirty'
SOURCE_PATHS = ("scripts", "src", "site", "report", "results", "docs", "README.md")
# environment variables that would point git at another repository, index or work tree
_LOCATION_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_COMMON_DIR", "GIT_OBJECT_DIRECTORY", "GIT_NAMESPACE")


def _env(**extra) -> dict:
    return {**{k: v for k, v in os.environ.items() if k not in _LOCATION_VARS}, **extra}


def git(where, *args, **env) -> str:
    """Run git in the directory ``where`` with extra environment variables; its output, RuntimeError on failure."""
    r = subprocess.run(["git", *args], cwd=where, env=_env(**env), stdin=subprocess.DEVNULL, capture_output=True,
                       encoding="utf-8", errors="replace")
    if r.returncode:
        raise RuntimeError(f"git {' '.join(args)}: {(r.stderr or r.stdout).strip()}")
    return r.stdout.rstrip()


def rev(repo, name: str) -> str:
    """The commit ``name`` points to, or '' if there is none."""
    try:
        return git(repo, "rev-parse", "--verify", "--quiet", f"{name}^{{commit}}")
    except RuntimeError:
        return ""


def tree_files(repo, treeish: str) -> list[str]:
    """Paths of every file in ``treeish``."""
    return [p for p in git(repo, "ls-tree", "-r", "-z", "--full-tree", "--name-only", treeish).split("\0") if p]


def source_commit(repo) -> tuple[str, list[str]]:
    """Short HEAD commit, with '+dirty' if SOURCE_PATHS differ from it (untracked files included), and the changes."""
    status = git(repo, "--no-optional-locks", "status", "--porcelain", "--untracked-files=all", "--", *SOURCE_PATHS)
    changes = [line for line in status.splitlines() if line]
    return git(repo, "rev-parse", "--short", "HEAD") + ("+dirty" if changes else ""), changes


def build(out: Path) -> str:
    """Build the site into the new directory ``out`` with scripts/build_site.py and return its output;
    RuntimeError if the build fails or has a broken local link or anchor."""
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "build_site.py"), "--out", str(out)], cwd=ROOT, env=_env(),
                       stdin=subprocess.DEVNULL, capture_output=True, text=True)
    log = (r.stdout + r.stderr).strip()
    if r.returncode:
        raise RuntimeError(f"the build failed (exit status {r.returncode}):\n{log}")
    problems = sb.check_links(out)
    if problems:
        raise RuntimeError(f"{len(problems)} broken links in the build:\n" + "\n".join(problems))
    if not (out / "index.html").is_file() or not (out / MARKER).is_file():
        raise RuntimeError(f"the build has no index.html or no {MARKER}: has build_site.py changed?")
    return log


def commit_pages(repo, build_dir, message: str, branch: str = BRANCH) -> str:
    """Commit every file in ``build_dir`` plus an empty .nojekyll as the next commit of ``branch``; returns it.

    A temporary index is used, so the working tree, the index, HEAD and the current branch of ``repo`` stay as
    they are. Refused if ``branch`` is checked out in a worktree or the commit would not hold exactly the files
    of ``build_dir``."""
    repo, build_dir = Path(repo).resolve(), Path(build_dir).resolve()
    ref = f"refs/heads/{branch}"
    if f"branch {ref}" in git(repo, "worktree", "list", "--porcelain").splitlines():
        raise RuntimeError(f"{branch} is checked out in a worktree; switch that worktree to another branch first")
    try:  # the repository's own identity (.git/config), never a global one
        name, email = git(repo, "config", "--local", "user.name"), git(repo, "config", "--local", "user.email")
    except RuntimeError:
        raise RuntimeError("set user.name and user.email in the repository's own config (git config --local ...)") from None
    (build_dir / ".nojekyll").write_bytes(b"")
    plumbing = ("--git-dir", git(repo, "rev-parse", "--absolute-git-dir"), "--work-tree", str(build_dir))
    with tempfile.TemporaryDirectory() as tmp:
        index = str(Path(tmp) / "index")
        # --force: the repository's ignore rules (e.g. /data in .git/info/exclude) must not drop files of the build
        git(build_dir, *plumbing, "add", "--all", "--force", "--", ".", GIT_INDEX_FILE=index)
        tree = git(build_dir, *plumbing, "write-tree", GIT_INDEX_FILE=index)
    on_disk = {unicodedata.normalize("NFC", p.relative_to(build_dir).as_posix()) for p in build_dir.rglob("*") if p.is_file()}
    in_tree = {unicodedata.normalize("NFC", p) for p in tree_files(repo, tree)}
    if on_disk != in_tree:
        raise RuntimeError(f"the commit would not hold exactly the build: left out {sorted(on_disk - in_tree)[:10]}, "
                           f"extra {sorted(in_tree - on_disk)[:10]}")
    old = rev(repo, ref)
    parent = old or rev(repo, f"refs/remotes/origin/{branch}")
    who = dict(GIT_AUTHOR_NAME=name, GIT_AUTHOR_EMAIL=email, GIT_COMMITTER_NAME=name, GIT_COMMITTER_EMAIL=email)
    commit = git(repo, "commit-tree", tree, *(["-p", parent] if parent else []), "-m", message, **who)
    git(repo, "update-ref", "-m", f"deploy_pages: {message.splitlines()[0]}", ref, commit, old)  # old '' = must not exist
    return commit


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--allow-dirty", action="store_true",
                    help=f"deploy although {', '.join(SOURCE_PATHS)} have uncommitted changes (the message says +dirty)")
    ap.add_argument("--push", action="store_true", help=f"then run `git push origin {BRANCH}` (a normal push, never forced)")
    args = ap.parse_args(argv)
    try:
        source, changes = source_commit(ROOT)
        if changes and not args.allow_dirty:
            raise RuntimeError(f"uncommitted changes in {', '.join(SOURCE_PATHS)} (commit them or pass --allow-dirty):\n"
                               + "\n".join(changes))
        message = (f"Pages: site built from {source} by scripts/build_site.py\n\nSource commit {git(ROOT, 'rev-parse', 'HEAD')}"
                   + (".\nUncommitted changes at build time:\n" + "\n".join(changes) if changes else "."))
        with tempfile.TemporaryDirectory(prefix="opmsquid-pages-") as tmp:
            out = Path(tmp) / "site"
            print(build(out))
            (out / MARKER).unlink()
            commit = commit_pages(ROOT, out, message)
    except RuntimeError as e:
        raise SystemExit(f"deploy_pages: {e}") from None
    print(f"{BRANCH} -> {commit} ({len(tree_files(ROOT, commit))} files; parent {rev(ROOT, commit + '^') or 'none'}): "
          f"{message.splitlines()[0]}")
    if not args.push:
        print(f"not pushed: run `git push origin {BRANCH}` (or deploy again with --push) to publish")
        return
    if subprocess.run(["git", "push", "origin", BRANCH], cwd=ROOT, env=_env()).returncode:
        raise SystemExit(f"deploy_pages: the push failed; the local {BRANCH} keeps {commit}")


if __name__ == "__main__":
    main()
