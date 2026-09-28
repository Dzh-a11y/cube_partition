"""Enumerate fixed-count cube partitions using one no-good cut per solution."""

import argparse
from dataclasses import asdict
import json
import math
import os
from pathlib import Path
import time

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, vstack

from solve import Cube, candidates, validate


def enumerate_partitions(n, count, output, time_limit=60, max_solutions=None,
                         resume=False, progress=False):
    if type(n) is not int or type(count) is not int or n < 1 or count < 1:
        raise ValueError("n and count must be positive integers")
    if not math.isfinite(time_limit) or time_limit <= 0:
        raise ValueError("time_limit must be positive and finite")
    if max_solutions is not None and (type(max_solutions) is not int or max_solutions < 1):
        raise ValueError("max_solutions must be a positive integer")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    header = {"record_type": "metadata", "format_version": 1,
              "grid_side": n, "cube_count": count, "mode": "all_feasible",
              "symmetry": "rotations and reflections count as distinct placements"}
    cubes = candidates(n, count)
    cube_index = {c: i for i, c in enumerate(cubes)}
    exclusions, seen = [], set()
    histograms = {}

    def check_solution(selected):
        histogram = validate(n, count, selected)
        key = tuple(sorted(cube_index[c] for c in selected))
        if key in seen:
            raise ValueError("Duplicate solution")
        seen.add(key)
        exclusions.append(key)
        histogram_key = ",".join(f"{s}:{num}" for s, num in histogram.items())
        histograms[histogram_key] = histograms.get(histogram_key, 0) + 1
        return histogram

    if resume:
        with output.open(encoding="utf-8") as stream:
            if json.loads(next(stream)) != header:
                raise ValueError("Checkpoint parameters or format do not match")
            for line in stream:
                record = json.loads(line)
                if record["record_type"] != "solution" or record["id"] != len(exclusions) + 1:
                    raise ValueError("Invalid checkpoint record")
                check_solution([Cube(**c) for c in record["cubes"]])
    else:
        # Never silently overwrite earlier solutions.
        with output.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(header) + "\n")
    resumed_count = len(exclusions)
    started = time.monotonic()
    status, message = "running", ""
    summary_path = output.with_suffix(".summary.json")

    def save_summary():
        report = {**header, "record_type": "summary", "scipy_version": scipy.__version__,
                  "status": status, "message": message, "complete": status == "exhausted",
                  "solutions_found": len(exclusions), "new_solutions": len(exclusions) - resumed_count,
                  "all_saved_solutions_verified": True,
                  "size_histograms": histograms, "candidate_count": len(cubes),
                  "elapsed_this_run_seconds": round(time.monotonic() - started, 3),
                  "solutions_file": str(output.resolve())}
        temporary = summary_path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(report, indent=2, ensure_ascii=False,
                                         allow_nan=False) + "\n", encoding="utf-8")
        temporary.replace(summary_path)
        return report

    save_summary()
    # All-feasible enumeration needs only x variables; y and the maximizing
    # objective remain in solve.py for finding one optimal partition.
    rows, cols = [], []
    for index, c in enumerate(cubes):
        cells = list(c.cells(n)) + [n ** 3]
        rows.extend(cells)
        cols.extend([index] * len(cells))
    base = coo_matrix((np.ones(len(rows)), (rows, cols)),
                      shape=(n ** 3 + 1, len(cubes))).tocsc()
    rhs = np.r_[np.ones(n ** 3), count]
    try:
        with output.open("a", encoding="utf-8") as stream:
            while True:
                if max_solutions is not None and len(exclusions) >= max_solutions:
                    status = "solution_limit"
                    break
                remaining = time_limit - (time.monotonic() - started)
                if remaining <= 0:
                    status = "time_limit"
                    break
                if not cubes:
                    status, message = "exhausted", "No legal candidates"
                    break
                matrix, lower, upper = base, rhs, rhs
                if exclusions:
                    cut_rows = np.repeat(np.arange(len(exclusions)), count)
                    cut_cols = np.array(exclusions, dtype=int).ravel()
                    cuts = coo_matrix((np.ones(len(cut_cols)), (cut_rows, cut_cols)),
                                      shape=(len(exclusions), len(cubes))).tocsc()
                    matrix = vstack([base, cuts], format="csc")
                    lower = np.r_[rhs, np.full(len(exclusions), -np.inf)]
                    upper = np.r_[rhs, np.full(len(exclusions), count - 1)]
                remaining = time_limit - (time.monotonic() - started)
                if remaining <= 0:
                    status = "time_limit"
                    break
                result = milp(np.zeros(len(cubes)), integrality=np.ones(len(cubes)),
                              bounds=Bounds(0, 1), constraints=LinearConstraint(matrix, lower, upper),
                              options={"time_limit": remaining, "mip_rel_gap": 0.0})
                message = result.message
                if result.status == 2:
                    status = "exhausted"
                    break
                if result.status not in (0, 1):
                    status = "solver_error"
                    break
                if result.x is not None:
                    if np.max(np.abs(result.x - np.rint(result.x))) > 1e-5:
                        raise ValueError("Nonintegral solver incumbent")
                    selected = [c for c, value in zip(cubes, result.x) if value > 0.5]
                    histogram = check_solution(selected)
                    record = {"record_type": "solution", "id": len(exclusions),
                              "distinct_sizes": len(histogram), "size_counts": histogram,
                              "cubes": [asdict(c) for c in selected]}
                    stream.write(json.dumps(record, separators=(",", ":")) + "\n")
                    stream.flush()
                    os.fsync(stream.fileno())
                    save_summary()
                    if progress and (len(exclusions) <= 5 or len(exclusions) % 10 == 0):
                        print(f"Verified and saved {len(exclusions)} solutions", flush=True)
                if result.status == 1:
                    status = "time_limit"
                    break
                if result.x is None:
                    status = "solver_error"
                    break
    except KeyboardInterrupt:
        status, message = "interrupted", "Saved solutions can be resumed"
    except Exception:
        status, message = "error", "Enumeration failed; inspect the exception before resuming"
        save_summary()
        raise
    return save_summary()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--count", type=int, default=49)
    parser.add_argument("--time-limit", type=float, default=60,
                        help="Time budget for this run in seconds (default: 60)")
    parser.add_argument("--max-solutions", type=int,
                        help="Optional total solution cap, including resumed solutions")
    parser.add_argument("--output", type=Path, required=True, help="Append-only JSONL checkpoint")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    try:
        report = enumerate_partitions(args.n, args.count, args.output, args.time_limit,
                                      args.max_solutions, args.resume, progress=True)
    except (ValueError, OSError) as error:
        parser.exit(2, f"{error}\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if report["status"] in ("solver_error", "error"):
        raise SystemExit(3)


if __name__ == "__main__":
    main()
