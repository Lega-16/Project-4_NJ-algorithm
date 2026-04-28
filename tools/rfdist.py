# rfdist.py <tree1> <tree2>
#
# Implementation of Day's algorihtm for computing the rf-distance between
# tree1 to tree2 over the same set of leaves, i.e. the number of splits not
# found in both trees.The two trees are read from the commandline and are
# assume to be in Newick-format. The Newick-parser from Biopyton
# (see https://biopython.org/wiki/Phylo) is used.
#
# Christian Storm Pedersen <cstorm@birc.au.dk>

import sys
from Bio import Phylo


def count_dublets(l):
    """
    returns the number of dublets in a list
    """
    return len(l) - len(list(set(l)))


def dfs(node, splits):
    """
    performs a dfs traversal of a Phylo tree from 'node' and adds splits to the
    list 'splits' that form a consecutive interval cf. the naming of the leaves
    specified in the dictionary 'dfs_num' that maps leaf names to numbers
    """
    if node.is_terminal():
        minval = maxval = dfs_num[node.name]
        size = 1
    else:
        minval = sys.maxsize
        maxval = -sys.maxsize
        size = 0
        for child in node.clades:
            child_min, child_max, child_size = dfs(child, splits)
            if child_min < minval:
                minval = child_min
            if child_max > maxval:
                maxval = child_max
            size = size + child_size
        if size == maxval - minval + 1:
            splits.append((minval, maxval))
    return minval, maxval, size


# Read the two trees from the commandline
tree1 = Phylo.read(sys.argv[1], "newick")
tree2 = Phylo.read(sys.argv[2], "newick")

# Reroot the two trees by 'outgrouping' the same leaf.
root_leaf = tree1.get_terminals()[0].name
tree1.root_with_outgroup(root_leaf)
tree2.root_with_outgroup(root_leaf)

# After rerooting the two trees with 'root_leaf' as outgroup, they look as:
#
#          tree.root
#      ________|_______
#     |               v
# root_leaf       ____|____
#                 |        |
#
# The number of 'non-trivial' splits is equal to the number of internal nodes
# in the subtree rooted at v, i.e. the number of internal nodes in the entire
# tree minus 2, since tree.root and v correspond to trivial splits in the
# original input tree.

# Make mapping from leaf names to dfs-numbers
dfs_num = {}
num = 1
for leaf in tree1.find_clades("", True, "postorder"):
    dfs_num[leaf.name] = num
    num = num + 1

# Collect all splits in tree1
splits = []
dfs(tree1.root, splits)

# Remove the two trivial splits that occur due to the rooting of tree1
splits = splits[:-2]

# Add potential shared splits form tree2
dfs(tree2.root, splits)

# Remove the two trivial splits that occur due to the rooting of tree2
splits = splits[:-2]

# The RF-distance is the number of unique splits in tree1 and tree2
num_of_nontrivial_splits = len(
    tree1.get_nonterminals()) - 2 + len(tree2.get_nonterminals()) - 2
num_of_shared_nontrivial_splits = count_dublets(splits)
rfdist = num_of_nontrivial_splits - 2 * num_of_shared_nontrivial_splits

print("RF Distance (number of unique splits):", rfdist)
print("Fraction of non-trivial splits:       ",
      rfdist / num_of_nontrivial_splits)
