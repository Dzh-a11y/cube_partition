import unittest

from solve import Cube, candidates, solve, validate


class CubePartitionTests(unittest.TestCase):
    def test_known_49_cube_partition(self):
        cubes = ([Cube(i, j, 0, 3) for i in (0, 3) for j in (0, 3)]
                 + [Cube(i, j, 3, 2) for i in (0, 2, 4) for j in (0, 2, 4)]
                 + [Cube(i, j, 5, 1) for i in range(6) for j in range(6)])
        self.assertEqual(validate(6, 49, cubes), {1: 36, 2: 9, 3: 4})

    def test_equal_volume_does_not_imply_valid_geometry(self):
        with self.assertRaisesRegex(ValueError, "Overlap"):
            validate(2, 8, [Cube(0, 0, 0, 1)] * 8)

    def test_outside_container(self):
        with self.assertRaisesRegex(ValueError, "outside"):
            validate(2, 1, [Cube(1, 0, 0, 2)])

    def test_candidate_boundaries(self):
        self.assertEqual(len(candidates(2, 1)), 9)
        self.assertTrue(all(max(c.i, c.j, c.k) + c.s <= 3 for c in candidates(3, 20)))

    def test_two_sizes_optimum(self):
        # 3^3 = one 2^3 cube + nineteen unit cubes; 20 cubes total.
        result = solve(3, 20, 10)
        self.assertTrue(result["optimal_on_grid"])
        self.assertTrue(result["verified"])
        self.assertEqual(result["distinct_sizes"], 2)
        self.assertEqual(result["size_counts"], {1: 19, 2: 1})

    def test_unit_cubes(self):
        result = solve(2, 8, 10)
        self.assertEqual(result["distinct_sizes"], 1)
        self.assertTrue(result["optimal_on_grid"])

    def test_infeasible_count(self):
        self.assertEqual(solve(2, 7, 10)["status"], "infeasible")
        self.assertEqual(solve(2, 9, 10)["status"], "infeasible")


if __name__ == "__main__":
    unittest.main()
