import subprocess
from pathlib import Path
import csv

BASE_DIR = Path(__file__).resolve().parent.parent

TREES_DIR = BASE_DIR / "results/trees"
RESULTS_FILE = BASE_DIR / "results" / "rf_distances.csv"

rows = []

files = sorted(TREES_DIR.glob("*_quicktree.newick"))

for qt_file in files:
    base = qt_file.name.replace("_quicktree.newick", "")

    qt = TREES_DIR / f"{base}_quicktree.newick"
    rnj = TREES_DIR / f"{base}_rapidnj.newick"
    our = TREES_DIR / f"{base}_our_nj.newick"

    def rf(a, b):
        print(f" Comparing {a.name} vs {b.name}")

        result = subprocess.check_output(
            ["python3", str(BASE_DIR / "tools" / "rfdist.py"), str(a), str(b)],
            text=True
        )

        for line in result.splitlines():
            if line.startswith("RF Distance"):
                return int(line.split(":")[1].strip())

        raise ValueError(f"Could not parse RF distance from output:\n{result}")

    rows.append({
        "matrix": base,
        "qt_vs_rnj": rf(qt, rnj),
        "qt_vs_our": rf(qt, our),
        "rnj_vs_our": rf(rnj, our),
    })

    print(f"Processing {base}")

with open(RESULTS_FILE, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["matrix", "qt_vs_rnj", "qt_vs_our", "rnj_vs_our"]
    )
    writer.writeheader()
    writer.writerows(rows)

print("RF distances saved in results/rf_distances.csv")
