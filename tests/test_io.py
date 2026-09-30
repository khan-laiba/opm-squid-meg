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

    def test_run_commit_format(self):
        self.assertRegex(io.RUN_COMMIT, r"^([0-9a-f]{7,40}(\+dirty)?|unknown)$")


if __name__ == "__main__":
    unittest.main()
