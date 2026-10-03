import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from opmsquid import io


class TestProvenance(unittest.TestCase):
    def test_run_commit_is_fixed_at_import(self):
        # outputs carry the commit the process loaded, not whatever HEAD is when they are written
        with tempfile.TemporaryDirectory() as d, mock.patch.object(io, "RUN_COMMIT", "abc1234"), \
                mock.patch.object(io, "git_commit", side_effect=AssertionError("must not be re-queried")):
            path = Path(d) / "out.json"
            io.write_json({"x": np.float64(1.5)}, path)
            prov = json.loads(path.read_text())["provenance"]
        self.assertEqual(prov["commit"], "abc1234")

    def test_csv_status_line_round_trip(self):
        import csv

        with tempfile.TemporaryDirectory() as d, mock.patch.object(io, "RUN_COMMIT", "abc1234"):
            path = Path(d) / "t.csv"
            with open(path, "w", newline="") as fh:
                io.csv_status(fh, "NEW (a test;\n two lines)")
                wr = csv.writer(fh)
                wr.writerow(["a", "b"])
                wr.writerow([1, "x,y"])
            first = path.read_text().splitlines()[0]
            rows = io.read_csv(path)
        self.assertEqual(first, "# NEW (a test; two lines) | commit abc1234")
        self.assertEqual(rows, [{"a": "1", "b": "x,y"}])

    def test_run_commit_format(self):
        self.assertRegex(io.RUN_COMMIT, r"^([0-9a-f]{7,40}(\+dirty)?|unknown)$")


if __name__ == "__main__":
    unittest.main()
