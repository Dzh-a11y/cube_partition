"""Exact-cover 0-1 MIP: partition a cube and maximize distinct side lengths."""

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


@dataclass(frozen=True)
class Cube:
    i: int
    j: int
    k: int
    s: int

    def cells(self, n):
        for x in range(self.i, self.i + self.s):
            for y in range(self.j, self.j + self.s):
                for z in range(self.k, self.k + self.s):
                    yield (x * n + y) * n + z


def candidates(n, count):
    # Every other selected cube has volume >= 1: this pruning is exact.
    return [
        Cube(i, j, k, s)
        for s in range(1, n + 1)
        if s ** 3 + count - 1 <= n ** 3
        for i in range(n - s + 1)
        for j in range(n - s + 1)
        for k in range(n - s + 1)
    ]


def validate(n, count, selected):
    """Independent integer geometry check, not solver residuals alone."""
    if len(selected) != count:
        raise ValueError("Wrong number of cubes")
    coverage = [0] * (n ** 3)
    for c in selected:
        values = (c.i, c.j, c.k, c.s)
        if any(type(v) is not int for v in values):
            raise ValueError("Coordinates and side lengths must be integers")
        if c.s < 1 or min(c.i, c.j, c.k) < 0:
            raise ValueError("Invalid cube")
        if max(c.i, c.j, c.k) + c.s > n:
            raise ValueError("Cube extends outside the container")
        for cell in c.cells(n):
            coverage[cell] += 1
    if any(value != 1 for value in coverage):
        raise ValueError("Overlap or uncovered cells")
    return dict(sorted(Counter(c.s for c in selected).items()))


def solve(n=6, count=49, time_limit=60.0):
    if type(n) is not int or type(count) is not int or n < 1 or count < 1:
        raise ValueError("n and count must be positive integers")
    if not math.isfinite(time_limit) or time_limit <= 0:
        raise ValueError("time_limit must be positive and finite")
    report = {"grid_side": n, "cube_count": count,
              "scipy_version": scipy.__version__,
              "scope": "Axis-aligned integer-grid partitions of [0,N]^3",
              "objective": "maximize number of distinct side lengths",
              "optimal_on_grid": False, "verified": False}
    if count > n ** 3:
        return {**report, "status": "infeasible", "message": "Too many cubes for this grid"}

    cubes = candidates(n, count)
    sizes = sorted({c.s for c in cubes})
    m, q, volume = len(cubes), len(sizes), n ** 3
    size_index = {s: t for t, s in enumerate(sizes)}
    rows, cols, values = [], [], []

    def add(row, col, value):
        rows.append(row)
        cols.append(col)
        values.append(value)

    # Rows: exact coverage; exact count; two linking inequalities per size.
    for index, c in enumerate(cubes):
        for cell in c.cells(n):
            add(cell, index, 1)
        add(volume, index, 1)
        t = size_index[c.s]
        add(volume + 1 + 2 * t, index, 1)
        add(volume + 2 + 2 * t, index, 1)
    lower = np.ones(volume + 1 + 2 * q)
    upper = np.ones_like(lower)
    lower[volume] = upper[volume] = count
    for t, s in enumerate(sizes):
        # sum(x of size s) - y_s >= 0
        row = volume + 1 + 2 * t
        add(row, m + t, -1)
        lower[row], upper[row] = 0, np.inf
        # sum(x of size s) - M_s*y_s <= 0
        add(row + 1, m + t, -min(count, volume // s ** 3))
        lower[row + 1], upper[row + 1] = -np.inf, 0
    matrix = coo_matrix((values, (rows, cols)),
                        shape=(len(lower), m + q)).tocsc()
    objective = np.r_[np.zeros(m), -np.ones(q)]  # SciPy minimizes.
    result = milp(objective, integrality=np.ones(m + q),
                  bounds=Bounds(0, 1),
                  constraints=LinearConstraint(matrix, lower, upper),
                  options={"time_limit": time_limit, "mip_rel_gap": 0.0})
    statuses = {0: "optimal", 1: "limit_reached", 2: "infeasible",
                3: "unbounded", 4: "solver_error"}
    report.update(status=statuses.get(result.status, "unknown"),
                  message=result.message, candidate_count=m,
                  binary_variable_count=m + q,
                  constraint_count=len(lower),
                  optimal_on_grid=result.status == 0)
    for key in ("mip_gap", "mip_node_count"):
        value = getattr(result, key, None)
        if value is not None and math.isfinite(float(value)):
            report[key] = float(value)
    dual = getattr(result, "mip_dual_bound", None)
    if dual is not None and math.isfinite(float(dual)):
        report["size_count_upper_bound"] = math.floor(-float(dual) + 1e-6)
    if result.x is not None:
        if np.max(np.abs(result.x - np.rint(result.x))) > 1e-5:
            raise ValueError("Solver returned a nonintegral incumbent")
        chosen = [c for c, value in zip(cubes, result.x[:m]) if value > 0.5]
        histogram = validate(n, count, chosen)
        active = {s for t, s in enumerate(sizes) if result.x[m + t] > 0.5}
        if active != set(histogram):
            raise ValueError("Size indicator mismatch")
        if abs(result.fun + len(histogram)) > 1e-5:
            raise ValueError("Objective mismatch")
        report.update(verified=True, distinct_sizes=len(histogram),
                      size_counts=histogram, cubes=[asdict(c) for c in chosen])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=6, help="Container grid side (default: 6)")
    parser.add_argument("--count", type=int, default=49)
    parser.add_argument("--time-limit", type=float, default=60000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.n < 1 or args.count < 1 or not math.isfinite(args.time_limit) or args.time_limit <= 0:
        parser.error("n, count and time-limit must be positive and finite")
    report = solve(args.n, args.count, args.time_limit)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2,
                                          allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "cubes"},
                     ensure_ascii=False, indent=2, allow_nan=False))
    if not report["verified"]:
        raise SystemExit(2 if report["status"] == "infeasible" else 3)


if __name__ == "__main__":
    main()
