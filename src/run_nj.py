
import sys
from nj_functions import run_nj

def main():
    if len(sys.argv) != 2:
        print("Usage: python run_nj.py <phylip_file>")
        sys.exit(1)

    infile = sys.argv[1]
    newick = run_nj(infile)
    print(newick)

if __name__ == "__main__":
    main()
