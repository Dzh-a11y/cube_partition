import json
from pathlib import Path
import tempfile
import unittest

from enumerate_solutions import enumerate_partitions


class EnumerationTests(unittest.TestCase):
    def test_all_eight_placements_and_resume(self):
        # A 2-cube has exactly 8 possible origins in a 3-cube. The remaining
        # nineteen cubes must be unit cubes, so this exhaustively checks cuts.
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "solutions.jsonl"
            first = enumerate_partitions(3, 20, output, 10, max_solutions=3)
            self.assertFalse(first["complete"])
            self.assertEqual(first["solutions_found"], 3)
            result = enumerate_partitions(3, 20, output, 10, resume=True)
            self.assertTrue(result["complete"])
            self.assertEqual(result["solutions_found"], 8)
            self.assertEqual(result["new_solutions"], 5)
            records = [json.loads(line) for line in output.read_text().splitlines()][1:]
            origins = {(c["i"], c["j"], c["k"]) for r in records for c in r["cubes"] if c["s"] == 2}
            self.assertEqual(origins, {(i, j, k) for i in range(2) for j in range(2) for k in range(2)})
            with self.assertRaises(FileExistsError):
                enumerate_partitions(3, 20, output, 10)
            with self.assertRaisesRegex(ValueError, "parameters"):
                enumerate_partitions(3, 19, output, 10, resume=True)

    def test_infeasible_and_single_solution(self):
        with tempfile.TemporaryDirectory() as folder:
            for count, expected in ((7, 0), (8, 1), (9, 0)):
                report = enumerate_partitions(2, count, Path(folder) / f"{count}.jsonl", 10)
                self.assertTrue(report["complete"])
                self.assertEqual(report["solutions_found"], expected)

    def test_time_limit_is_not_exhaustion(self):
        with tempfile.TemporaryDirectory() as folder:
            report = enumerate_partitions(3, 20, Path(folder) / "limited.jsonl", 1e-12)
            self.assertEqual(report["status"], "time_limit")
            self.assertFalse(report["complete"])


if __name__ == "__main__":
    unittest.main()
