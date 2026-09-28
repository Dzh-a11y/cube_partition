"""Select diverse saved placements, treating all 48 cube symmetries as equivalent."""
import itertools
import json
from pathlib import Path

from solve import Cube, validate

ROOT = Path(__file__).resolve().parents[2]
source = ROOT / "outputs/cube_partition_mip/feasible_n6.jsonl"
records = [json.loads(line) for line in source.read_text().splitlines()]
assert records[0]["grid_side"] == 6 and records[0]["cube_count"] == 49
transforms = list(itertools.product(itertools.permutations(range(3)),
                                    itertools.product((0, 1), repeat=3)))


def orbit(cubes):
    return [frozenset(tuple(6 - c[p[a]] - c[3] if flips[a] else c[p[a]]
                           for a in range(3)) + (c[3],) for c in cubes)
            for p, flips in transforms]


unique = {}
for record in records[1:]:
    validate(6, 49, [Cube(**c) for c in record["cubes"]])
    cubes = [tuple(c[k] for k in ("i", "j", "k", "s")) for c in record["cubes"]]
    images = orbit(cubes)
    key = min(tuple(sorted(image)) for image in images)
    if key not in unique:
        unique[key] = {"id": record["id"], "cubes": cubes, "orbit": images}

pool = list(unique.values())
chosen = [pool.pop(0)]
while len(chosen) < 10:
    # Maximize the minimum number of changed candidate cubes, allowing the
    # best alignment under rotations AND reflections for each comparison.
    def distance(item):
        return min(49 - max(len(set(other["cubes"]) & image) for image in item["orbit"])
                   for other in chosen)
    next_item = max(pool, key=lambda item: (distance(item), -item["id"]))
    pool.remove(next_item)
    chosen.append(next_item)

payload = {"n": 6, "count": 49, "source_solution_count": len(records) - 1,
           "symmetry_classes": len(unique), "selection": "greedy maximin under all 48 cube symmetries",
           "solutions": [{"id": c["id"], "cubes": c["cubes"]} for c in chosen]}
destination = ROOT / "outputs/cube_partition_mip/demo_selected.json"
destination.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps({"selected_ids": [x["id"] for x in chosen],
                  "symmetry_classes": len(unique), "output": str(destination)}))
