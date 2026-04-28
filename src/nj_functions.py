import math


def read_phylip_distance_matrix(filepath):
    with open(filepath, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    n = int(lines[0])
    names = []
    D = []

    for line in lines[1:1 + n]:
        parts = line.split()
        names.append(parts[0])
        D.append([float(x) for x in parts[1:]])

    return names, D


def neighbor_joining(names, D):
    current_names = names[:]
    current_D = [row[:] for row in D]

    while len(current_names) > 2:
        n = len(current_names)

        r = [sum(current_D[i]) / (n - 2) for i in range(n)]

        min_value = math.inf
        min_i = None
        min_j = None

        for i in range(n):
            for j in range(i + 1, n):
                value = current_D[i][j] - r[i] - r[j]

                if value < min_value:
                    min_value = value
                    min_i = i
                    min_j = j

        i = min_i
        j = min_j
        dij = current_D[i][j]

        length_i = 0.5 * (dij + r[i] - r[j])
        length_j = dij - length_i

        new_node = (
            f"({current_names[i]}:{length_i:.10f},"
            f"{current_names[j]}:{length_j:.10f})"
        )

        keep = [k for k in range(n) if k not in (i, j)]

        new_D = []
        new_names = []

        for a in keep:
            row = []
            for b in keep:
                row.append(current_D[a][b])
            new_D.append(row)
            new_names.append(current_names[a])

        new_distances = []

        for k in keep:
            d_new_k = 0.5 * (
                current_D[i][k] + current_D[j][k] - dij
            )
            new_distances.append(d_new_k)

        for row, distance in zip(new_D, new_distances):
            row.append(distance)

        new_distances.append(0.0)
        new_D.append(new_distances)
        new_names.append(new_node)

        current_D = new_D
        current_names = new_names

    final_length = current_D[0][1] / 2

    newick = (
        f"({current_names[0]}:{final_length:.10f},"
        f"{current_names[1]}:{final_length:.10f});"
    )

    return newick


def run_nj(filepath):
    names, D = read_phylip_distance_matrix(filepath)
    return neighbor_joining(names, D)
