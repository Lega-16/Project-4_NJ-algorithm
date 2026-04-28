import sys
from nj_functions import run_nj


def main():
    if len(sys.argv) != 3:
        print("Usage: python src/run_nj.py input.phy output.newick")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    newick = run_nj(input_file)

    with open(output_file, "w") as f:
        f.write(newick + "\n")


if __name__ == "__main__":
    main()
