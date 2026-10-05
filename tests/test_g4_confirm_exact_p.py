"""The exact sign-flip p of scripts/g4_confirm_exact_p.py: against brute-force enumeration and against the enumeration
that opmsquid.detection.sign_flip_p performs itself up to 20 non-zero differences."""
import importlib.util
import itertools
import random
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from opmsquid import detection

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("g4_confirm_exact_p", REPO / "scripts" / "g4_confirm_exact_p.py")
E = importlib.util.module_from_spec(spec)
spec.loader.exec_module(E)


class TestExactSignFlip(unittest.TestCase):
    def test_brute_force(self):
        rng = random.Random(7)
        for _ in range(200):
            d = [rng.choice([-4, -2, -1, 0, 1, 2, 3, 5]) for _ in range(rng.randint(1, 11))]
            x = [v for v in d if v]
            p, n = E.exact_sign_flip(d)
            self.assertEqual(n, len(x))
            if not x:
                self.assertEqual(p, 1)
                continue
            hits = sum(abs(sum(s * v for s, v in zip(signs, x))) >= abs(sum(x))
                       for signs in itertools.product((-1, 1), repeat=len(x)))
            self.assertEqual(p, Fraction(hits, 2 ** len(x)), d)

    def test_matches_the_runs_own_enumeration(self):
        rng = np.random.default_rng(3)
        for _ in range(20):
            d = rng.integers(-6, 7, size=rng.integers(5, 21)).tolist()
            p, n = E.exact_sign_flip(d)
            if n <= 20:
                self.assertAlmostEqual(float(p), detection.sign_flip_p(d), places=12)

    def test_all_one_sided(self):
        p, n = E.exact_sign_flip([1] * 30)  # every pattern but the two extremes is less extreme
        self.assertEqual((p, n), (Fraction(2, 2 ** 30), 30))

    def test_holm(self):
        self.assertEqual(E.holm({"a": 0.01, "b": 0.04, "c": 0.03}), {"a": 0.03, "b": 0.06, "c": 0.06})


if __name__ == "__main__":
    unittest.main()
