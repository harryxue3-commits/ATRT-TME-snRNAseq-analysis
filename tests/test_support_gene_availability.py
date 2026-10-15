"""Regression check: an absent gene must not dilute an available-gene mean."""
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "03_astrocyte_support_analysis/analyze_astrocyte_support_modules.py"
SPEC = importlib.util.spec_from_file_location("support", SCRIPT)
support = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(support)


class AvailableGeneMeanTest(unittest.TestCase):
    def test_absent_gene_is_excluded_but_measured_zero_is_retained(self):
        with tempfile.TemporaryDirectory() as directory:
            matrix = Path(directory) / "counts.tsv.gz"
            with gzip.open(matrix, "wt") as handle:
                handle.write("c1\tc2\nG1\t4\t0\nG1\t1\t0\nG2\t0\t0\nOTHER\t5\t10\n")
            counts, libraries, _, duplicates = support.aggregate_counts(
                matrix, {"c1": "A", "c2": "T"}, {"G1", "G2", "ABSENT"}
            )
        self.assertEqual(set(counts.index), {"G1", "G2"})
        self.assertEqual(duplicates, ["G1"])
        self.assertEqual(counts.loc["G1", "A"], 5)
        self.assertEqual(libraries.to_dict(), {"A": 10, "T": 10})
        score = np.log2(counts.div(libraries, axis=1) * 1e6 + 1).mean(axis=0)
        self.assertAlmostEqual(score["A"], np.log2(500001) / 2)
        self.assertEqual(score["T"], 0)


if __name__ == "__main__":
    unittest.main()
