from pathlib import Path
from fractions import Fraction
import csv
import json
import sys
import tempfile
import unittest
import numpy as np
root = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(root / "analysis"), str(root / "src")]
from simulate import ceil_ticks, observe, generate
from bootstrap import Margin, daly_order, bootstrap
from storage import Store
from encoding import clean
from estimation import evaluate

class AnalysisExamples(unittest.TestCase):
    def test_represented_boundaries(self):
        count = 0
        for den in (4, 10):
            for k in range(41):
                self.assertEqual(int(ceil_ticks(Fraction(k, den), den)), k)
                t = float(Fraction(k, den))
                values = [t, np.nextafter(t, np.inf)]
                if k:
                    values.append(np.nextafter(t, -np.inf))
                for v in values:
                    a, b = float(v).as_integer_ratio()
                    expected = -((-a * den) // b)
                    self.assertEqual(int(ceil_ticks(v, den)), expected)
                    self.assertEqual(int(ceil_ticks(np.array([v]), den)[0]), expected)
                    count += 1
        self.assertEqual(count, 244)

    def test_event_before_censor(self):
        p, d = observe(np.array([1, 2]), np.array([1, 3]), np.array([1, 1]), np.array([1, 4]))
        np.testing.assert_array_equal(p, [[1, 1], [1, 0]])
        np.testing.assert_array_equal(d, [[1, 1], [3, 1]])

    def test_existing_index_paths(self):
        p = [(1, 1), (2, 1), (3, 1), (3, 0), (4, 1), (6, 0)]
        d = [(4, 1), (5, 1), (6, 0), (3, 0), (6, 0), (6, 0)]
        with (root / "examples/bootstrap_index_paths.csv").open() as f:
            for row in csv.DictReader(f):
                ip = np.array(list(map(int, row["p_indices"].split(";"))))
                idd = np.array(list(map(int, row["d_indices"].split(";"))))
                a = Margin(p, 6).curve(np.bincount(ip, minlength=6))
                b = Margin(d, 6).curve(np.bincount(idd, minlength=6))
                self.assertLess(abs(a[0] - float(row["rmst_p"])), 1e-10)
                self.assertLess(abs(b[0] - float(row["rmst_d"])), 1e-10)
                self.assertEqual(daly_order(*a[1:3], *b[1:3]), row["accepted"] == "TRUE")

    def test_empty_death_margin_rejection_is_not_an_order_violation(self):
        p = [(1, 1), (2, 1), (4, 0)]
        d = [(4, 0), (4, 0), (4, 0)]
        mp, tp, sp, sup_p = Margin(p, 3).curve(np.ones(3))
        md, td, sd, sup_d = Margin(d, 3).curve(np.ones(3))
        self.assertTrue(sup_p and sup_d)
        self.assertEqual(len(td), 0)
        self.assertTrue(np.all(sp <= 1))  # The death curve is identically one.
        self.assertFalse(daly_order(tp, sp, td, sd))

    def test_caps_and_support_retained(self):
        p = [(3, 0), (3, 0)]
        indices = (np.array([0, 1]), np.array([0, 1]))
        b = bootstrap(p, p, 2, 1, None, "literal", 2, 3, iter([indices] * 3))
        self.assertEqual((b["status"], b["attempted"], b["accepted"]), ("CAPPED", 3, 0))
        b = bootstrap(p, p, 2, 1, None, "paired", 2, 3, iter([indices] * 3))
        self.assertTrue(b["complete"])
        self.assertIsNone(b["correlation"])
        self.assertEqual(b["variance_difference"], 0)

    def test_exact_checkpoint(self):
        with tempfile.TemporaryDirectory() as folder:
            s = Store(Path(folder) / "examples.sqlite", "example", [(0, 0), (0, 1)])
            s.start(0, 0)
            s.finish(0, 0, "COMPLETED", clean({"exact": Fraction(1, 3)}))
            self.assertEqual(s.get(0, 0)[1]["exact"], "1/3")
            with self.assertRaises(ValueError):
                s.start(0, 0)
            s.start(0, 1)
            s.finish(0, 1, "CALCULATION_FAILED", {"error": "example"})
            self.assertEqual(s.get(0, 1)[0], "CALCULATION_FAILED")
            s.close()

    def test_joint_comparators(self):
        p = [(1, 1), (1, 1), (3, 0)]
        d = [(1, 1), (2, 1), (3, 0)]
        x = evaluate(p, d, 2, check_optimality=True)
        self.assertLessEqual(x["unrestricted_lower"], x["oracle"])
        self.assertLessEqual(x["oracle"], x["unrestricted_upper"])
        for release in ("separate_R1", "separate_R2"):
            b = x["levels"][release]
            self.assertLessEqual(b["lower"], x["oracle"] + 1e-12)
            self.assertGreaterEqual(b["upper"], x["oracle"] - 1e-12)

if __name__ == "__main__":
    unittest.main(verbosity=2)
