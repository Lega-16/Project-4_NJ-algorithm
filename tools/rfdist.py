import sys
from pathlib import Path


def parse_newick(text):
    text = text.strip().rstrip(";")
    i = 0
    graph = {}
    leaves = set()
    node_id = 0

    def new_internal():
        nonlocal node_id
        name = f"__internal_{node_id}"
        node_id += 1
        graph.setdefault(name, [])
        return name

    def add_edge(a, b):
        graph.setdefault(a, []).append(b)
        graph.setdefault(b, []).append(a)

    def skip_spaces():
        nonlocal i
        while i < len(text) and text[i].isspace():
            i += 1

    def skip_branch_length():
        nonlocal i
        skip_spaces()
        if i < len(text) and text[i] == ":":
            i += 1
            while i < len(text) and text[i] not in ",();":
                i += 1

    def parse_label():
        nonlocal i
        skip_spaces()

        if text[i] in "'\"":
            quote = text[i]
            i += 1
            start = i
            while i < len(text) and text[i] != quote:
                i += 1
            label = text[start:i]
            i += 1
        else:
            start = i
            while i < len(text) and text[i] not in ":,();":
                i += 1
            label = text[start:i].strip()

        label = label.strip().strip("'\"")
        skip_branch_length()

        graph.setdefault(label, [])
        leaves.add(label)
        return label

    def parse_subtree():
        nonlocal i
        skip_spaces()

        if text[i] == "(":
            node = new_internal()
            i += 1

            while True:
                child = parse_subtree()
                add_edge(node, child)

                skip_spaces()

                if i < len(text) and text[i] == ",":
                    i += 1
                    continue

                if i < len(text) and text[i] == ")":
                    i += 1
                    break

            skip_spaces()

            # ignore possible internal node label
            while i < len(text) and text[i] not in ":,();":
                i += 1

            skip_branch_length()
            return node

        return parse_label()

    parse_subtree()
    return graph, leaves


def leaves_on_side(graph, leaves, start, blocked):
    stack = [start]
    seen = set()
    found = set()

    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)

        if node in leaves:
            found.add(node)

        for nxt in graph[node]:
            if (node, nxt) == blocked or (nxt, node) == blocked:
                continue
            stack.append(nxt)

    return found


def splits_from_tree(path):
    graph, leaves = parse_newick(Path(path).read_text())
    all_leaves = set(leaves)
    n = len(all_leaves)
    splits = set()

    for u in graph:
        for v in graph[u]:
            if str(u) > str(v):
                continue

            side = leaves_on_side(graph, all_leaves, u, (u, v))

            if 1 < len(side) < n - 1:
                other = all_leaves - side

                if len(other) < len(side):
                    side = other
                elif len(other) == len(side) and sorted(other) < sorted(side):
                    side = other

                splits.add(frozenset(side))

    return splits, all_leaves


def rf_distance(file1, file2):
    splits1, leaves1 = splits_from_tree(file1)
    splits2, leaves2 = splits_from_tree(file2)

    if leaves1 != leaves2:
        raise ValueError("The two trees do not contain the same leaves")

    return len(splits1.symmetric_difference(splits2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 rfdist.py tree1.newick tree2.newick")
        sys.exit(1)

    print(rf_distance(sys.argv[1], sys.argv[2]))
