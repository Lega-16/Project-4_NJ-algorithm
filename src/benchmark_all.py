import csv
import subprocess
import time
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

DATA_DIR = BASE_DIR / "data" / "unique_distance_matrices"
RESULTS_DIR = BASE_DIR / "results"
TREES_DIR = BASE_DIR / "results" / "trees"
TOOLS_DIR = BASE_DIR / "tools"

QUICKTREE = BASE_DIR / "quicktree-2.5" / "quicktree-2.5" / "quicktree"
RAPIDNJ = BASE_DIR / "rapidNJ" / "rapidNJ-master" / "bin" / "rapidnj"
OUR_NJ_SCRIPT = BASE_DIR / "src" / "run_nj.py"

RESULTS_DIR.mkdir(exist_ok=True)
TREES_DIR.mkdir(exist_ok=True)


def run_and_time(command, output_path):
    start = time.perf_counter()

    with open(output_path, "w") as out:
        subprocess.run(
            command,
            stdout=out,
            stderr=subprocess.DEVNULL,
            check=True
        )

    end = time.perf_counter()
    return end - start


rows = []

phy_files = sorted(DATA_DIR.glob("*.phy"))

for phy_file in phy_files:
    matrix_name = phy_file.stem
    print(f"Processing {phy_file.name}")

    quicktree_tree = TREES_DIR / f"{matrix_name}_quicktree.newick"
    rapidnj_tree = TREES_DIR / f"{matrix_name}_rapidnj.newick"
    our_tree = TREES_DIR / f"{matrix_name}_our_nj.newick"

    quicktree_time = run_and_time(
        [str(QUICKTREE), "-in", "m", str(phy_file)],
        quicktree_tree
    )

    rapidnj_time = run_and_time(
        [str(RAPIDNJ), str(phy_file)],
        rapidnj_tree
    )

    our_nj_time = run_and_time(
        ["python3", str(OUR_NJ_SCRIPT), str(phy_file), str(our_tree)],
        our_tree
    )

    rows.append({
        "matrix": phy_file.name,
        "quicktree_time": quicktree_time,
        "rapidnj_time": rapidnj_time,
        "our_nj_time": our_nj_time,
        "speedup_vs_quicktree": quicktree_time / our_nj_time,
        "speedup_vs_rapidnj": rapidnj_time / our_nj_time
    })

output_csv = RESULTS_DIR / "timings.csv"

with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "matrix",
            "quicktree_time",
            "rapidnj_time",
            "our_nj_time",
            "speedup_vs_quicktree",
            "speedup_vs_rapidnj"
        ]
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Done. Results saved in {output_csv}")
