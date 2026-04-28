from pathlib import Path
import re


def leaves(path):
    text = Path(path).read_text()
    text = re.sub(r":[^,();]+", "", text)
    tokens = re.findall(r"[A-Za-z0-9_.|/-]+", text)
    return set(tokens)


qt = leaves("results/trees/1347_FAINT_quicktree.newick")
rnj = leaves("results/trees/1347_FAINT_rapidnj.newick")

print("QuickTree leaves:", len(qt))
print("RapidNJ leaves:", len(rnj))

print("Only in QuickTree:")
print(sorted(qt - rnj)[:30])

print("Only in RapidNJ:")
print(sorted(rnj - qt)[:30])
