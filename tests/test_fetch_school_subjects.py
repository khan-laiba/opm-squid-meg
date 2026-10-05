"""scripts/fetch_school_subjects.py, offline: both manifests (the surfaces and the MRI quality check) read into one form
with the exact S3 object version of every file, malformed entries refused, the size and SHA-256 check on temporary
files, a download that differs from its entry or breaks off discarded, and --check, which downloads nothing. No test
reaches the network: urllib.request.urlopen is replaced by a stub that serves bytes from memory."""
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import tempfile
import tomllib
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
FREESURFER = "https://s3.amazonaws.com/openneuro.org/ds005234/derivatives/freesurfer/"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


F = _load("fetch_school_subjects")


def raw_entry(dest, content, version="v1", qc=False):
    """A manifest entry for `content`, in the form of the quality-check manifest (s3:// URI, version_id) or of the
    surfaces' manifest (https address, version)."""
    e = {"dest": dest, "size": len(content), "sha256": hashlib.sha256(content).hexdigest()}
    if qc:
        return {**e, "source": f"s3://example-bucket/{dest}", "version_id": version}
    return {**e, "source": f"https://s3.amazonaws.com/example-bucket/{dest}", "version": version}


class TestProjectManifests(unittest.TestCase):
    """configs/school_subjects_manifest.json and configs/school_subjects_qc_manifest.json as the fetch reads them."""

    @classmethod
    def setUpClass(cls):
        cls.surf = F.load(F.MANIFEST)
        cls.qc = F.load(F.QC_MANIFEST)

    def test_counts(self):
        self.assertEqual([k for k, _ in F.MANIFESTS], ["surfaces", "qc"])
        self.assertEqual((len(self.surf["files"]), sum(e["size"] for e in self.surf["files"])), (45, 124_961_312))
        self.assertEqual((len(self.qc["files"]), sum(e["size"] for e in self.qc["files"])), (6, 9_855_337))

    def test_exact_object_versions(self):
        for path, m in ((F.MANIFEST, self.surf), (F.QC_MANIFEST, self.qc)):
            raw = json.loads(path.read_text())["files"]
            for e, r in zip(m["files"], raw, strict=True):
                with self.subTest(manifest=path.name, dest=e["dest"]):
                    source = r["source"].replace("s3://", "https://s3.amazonaws.com/", 1)
                    self.assertTrue(source.startswith(FREESURFER))
                    self.assertEqual(e["url"], f"{source}?versionId={r.get('version', r.get('version_id'))}")
                    self.assertEqual((e["dest"], e["size"], e["sha256"]), (r["dest"], r["size"], r["sha256"]))

    def test_qc_covers_the_children(self):
        children = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())["anatomy"]["school"]
        self.assertEqual({e["dest"] for e in self.qc["files"]},
                         {f"{c['subject']}/mri/{name}" for c in children for name in ("T1.mgz", "seghead.mgz")})

    def test_no_destination_twice(self):
        dests = [e["dest"] for m in (self.surf, self.qc) for e in m["files"]]
        self.assertEqual(len(set(dests)), len(dests))


class TestMalformedManifests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = Path(tmp.name) / "manifest.json"

    def load(self, manifest):
        self.path.write_text(json.dumps(manifest))
        return F.load(self.path)

    def test_entries_refused(self):
        good = raw_entry("sub-X/mri/T1.mgz", b"abc", qc=True)
        self.assertEqual(self.load({"files": [good]})["files"][0]["url"],
                         "https://s3.amazonaws.com/example-bucket/sub-X/mri/T1.mgz?versionId=v1")
        cases = {
            "absolute destination": {**good, "dest": "/tmp/T1.mgz"},
            "destination outside the folder": {**good, "dest": "../T1.mgz"},
            "destination leaving through a subfolder": {**good, "dest": "sub-X/../../T1.mgz"},
            "empty destination": {**good, "dest": ""},
            "no size": {k: v for k, v in good.items() if k != "size"},
            "negative size": {**good, "size": -1},
            "size as text": {**good, "size": "3"},
            "short checksum": {**good, "sha256": good["sha256"][:63]},
            "upper-case checksum": {**good, "sha256": good["sha256"].upper()},
            "no object version": {k: v for k, v in good.items() if k != "version_id"},
            "empty object version": {**good, "version_id": ""},
            "plain http": {**good, "source": "http://s3.amazonaws.com/example-bucket/T1.mgz"},
            "local file": {**good, "source": "file:///etc/hosts"},
            "address with a query": {**good, "source": "https://s3.amazonaws.com/example-bucket/T1.mgz?versionId=x"},
        }
        for name, e in cases.items():
            with self.subTest(name), self.assertRaises(ValueError):
                self.load({"files": [e]})

    def test_totals_and_duplicates(self):
        a, b = raw_entry("a", b"aa"), raw_entry("b", b"bbb")
        self.assertEqual(len(self.load({"n_files": 2, "total_bytes": 5, "files": [a, b]})["files"]), 2)
        cases = {"file count": {"n_files": 3, "files": [a, b]},
                 "byte total": {"total_bytes": 6, "files": [a, b]},
                 "destination twice": {"files": [a, {**b, "dest": "a"}]},
                 "no file list": {"n_files": 0}}
        for name, m in cases.items():
            with self.subTest(name), self.assertRaises(ValueError):
                self.load(m)


class TestIntact(unittest.TestCase):
    def test_size_and_checksum(self):
        content = os.urandom(2_500_000)  # read in three of the 1-MiB blocks
        e = F.entry(raw_entry("T1.mgz", content))
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "T1.mgz"
            self.assertFalse(F.intact(p, e))  # missing
            p.write_bytes(content)
            self.assertEqual(F.sha256(p), hashlib.sha256(content).hexdigest())
            self.assertTrue(F.intact(p, e))
            changed = bytearray(content)
            changed[1_234_567] ^= 1
            p.write_bytes(bytes(changed))  # same size, one bit changed
            self.assertFalse(F.intact(p, e))
            p.write_bytes(content + b"\0")
            self.assertFalse(F.intact(p, e))
            p.unlink()
            p.mkdir()  # a folder in its place
            self.assertFalse(F.intact(p, e))


class _BreaksOff(io.BytesIO):
    """A response body that breaks off after its first bytes."""

    def read(self, n=-1):
        if self.tell():
            raise ConnectionResetError("connection reset (test)")
        return super().read(10)


class _Served:
    """Stands in for urllib.request.urlopen: serves bytes by URL from memory and records every request; a request for
    any other URL fails the test."""

    def __init__(self):
        self.body, self.requests, self.broken = {}, [], set()

    def __call__(self, url, timeout=None):
        self.requests.append(url)
        if url not in self.body:
            raise AssertionError(f"unexpected request: {url}")
        return (_BreaksOff if url in self.broken else io.BytesIO)(self.body[url])


class TestFetchAndCheck(unittest.TestCase):
    """main() on two temporary manifests (one in each form) and a temporary data folder."""

    FILES = {"surfaces": {"sub-A/surf/lh.white": b"white" * 3000, "sub-A/label/lh.aparc.annot": b"annot" * 10},
             "qc": {"sub-A/mri/T1.mgz": b"T1" * 5000, "sub-A/mri/seghead.mgz": b"mask" * 20}}

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.root = self.tmp / "external" / "school_subjects"
        self.served, self.urls, manifests = _Served(), {}, []
        for key, files in self.FILES.items():
            entries = []
            for i, (dest, content) in enumerate(files.items()):
                entries.append(raw_entry(dest, content, version=f"{key}.{i}", qc=key == "qc"))
                self.urls[dest] = f"https://s3.amazonaws.com/example-bucket/{dest}?versionId={key}.{i}"
                self.served.body[self.urls[dest]] = content
            path = self.tmp / f"{key}.json"
            path.write_text(json.dumps({"dataset": "Test dataset", "snapshot": "1.0", "citation": "Test 2026",
                                        "license": "CC0", "files": entries}))
            manifests.append((key, path))
        for patch in (mock.patch.object(F, "MANIFESTS", tuple(manifests)),
                      mock.patch.object(F.paths, "EXTERNAL", self.tmp / "external"),
                      mock.patch.object(F.urllib.request, "urlopen", self.served)):
            patch.start()
            self.addCleanup(patch.stop)

    def run_main(self, *argv):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            status = F.main(list(argv))
        return status, out.getvalue()

    def assert_no_part_files(self):
        self.assertEqual(list(self.tmp.rglob("*.part")), [])

    def test_check_downloads_nothing(self):
        status, out = self.run_main("--check")
        self.assertEqual(status, 1)
        self.assertEqual(out.count("0 of 2 files present and intact"), 2)
        for dest in self.urls:
            self.assertIn(f"missing or changed: {dest}", out)
        self.assertEqual(self.served.requests, [])
        self.assertFalse(self.root.exists())

    def test_fetches_what_is_missing_or_changed(self):
        kept = self.root / "sub-A/surf/lh.white"
        kept.parent.mkdir(parents=True)
        kept.write_bytes(self.FILES["surfaces"]["sub-A/surf/lh.white"])
        changed = self.root / "sub-A/mri/T1.mgz"
        changed.parent.mkdir(parents=True)
        changed.write_bytes(b"T2" * 5000)  # same size, other content: fetched again
        status, out = self.run_main()
        self.assertEqual(status, 0)
        self.assertIn("downloading 3 files", out)
        fetched = ("sub-A/label/lh.aparc.annot", "sub-A/mri/T1.mgz", "sub-A/mri/seghead.mgz")
        self.assertEqual(self.served.requests, [self.urls[d] for d in fetched])
        for files in self.FILES.values():
            for dest, content in files.items():
                self.assertEqual((self.root / dest).read_bytes(), content)
        self.assert_no_part_files()
        self.assertEqual(self.run_main("--check")[0], 0)
        self.assertEqual(self.run_main()[0], 0)  # nothing left to fetch
        self.assertEqual(len(self.served.requests), 3)

    def test_without_the_quality_check(self):
        self.assertEqual(self.run_main("--no-qc")[0], 0)
        self.assertEqual(self.served.requests, [self.urls[d] for d in self.FILES["surfaces"]])
        self.assertFalse((self.root / "sub-A/mri").exists())
        self.assertEqual(self.run_main("--check", "--no-qc")[0], 0)
        status, out = self.run_main("--check")
        self.assertEqual(status, 1)
        self.assertIn("2 of 2 files present and intact", out)
        self.assertIn("0 of 2 files present and intact", out)

    def test_download_that_differs_is_deleted(self):
        self.served.body[self.urls["sub-A/mri/T1.mgz"]] = b"T2" * 5000  # another object than the listed one
        with self.assertRaisesRegex(RuntimeError, "sub-A/mri/T1.mgz"):
            self.run_main()
        self.assertFalse((self.root / "sub-A/mri/T1.mgz").exists())
        self.assert_no_part_files()

    def test_download_that_breaks_off_leaves_nothing(self):
        self.served.broken.add(self.urls["sub-A/surf/lh.white"])
        with self.assertRaises(ConnectionResetError):
            self.run_main()
        self.assertFalse((self.root / "sub-A/surf/lh.white").exists())
        self.assert_no_part_files()

    def test_destination_in_both_manifests(self):
        qc = dict(F.MANIFESTS)["qc"]
        m = json.loads(qc.read_text())
        m["files"][0]["dest"] = "sub-A/surf/lh.white"
        qc.write_text(json.dumps(m))
        with self.assertRaises(ValueError):
            self.run_main("--check")
        self.assertEqual(self.run_main("--check", "--no-qc")[0], 1)  # the surfaces alone: listed once, not fetched yet
        self.assertEqual(self.served.requests, [])


if __name__ == "__main__":
    unittest.main()
