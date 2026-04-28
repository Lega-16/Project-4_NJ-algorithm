
# --- TUS FUNCIONES AQUÍ ---
# Pega aquí:



# compute_Q_matrix
# find_min_Q
# update_distance_matrix
# build_newick
# etc.

def read_phylip_distance_matrix(filepath):
    with open(filepath, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # Número de taxones
    n = int(lines[0])

    names = []
    matrix = []

    for line in lines[1:1+n]:
        parts = line.split()
        name = parts[0]
        distances = list(map(float, parts[1:]))
        names.append(name)
        matrix.append(distances)

    return names, matrix


def compute_r_values(D):
    n = len(D)
    r = []
    for i in range(n):
        r_i = sum(D[i]) / (n - 2)
        r.append(r_i)
    return r


import math

def compute_N_matrix(D, r):
    n = len(D)
    N = [[None]*n for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            if i == j:
                N[i][j] = math.inf
            else:
                N[i][j] = D[i][j] - (r[i] + r[j])
    return N


def find_min_N(N):
    n = len(N)
    min_val = math.inf
    min_i, min_j = None, None
    
    for i in range(n):
        for j in range(i+1, n):
            if N[i][j] < min_val:
                min_val = N[i][j]
                min_i, min_j = i, j
    
    return min_i, min_j


def compute_branch_lengths(D, r, i, j):
    dij = D[i][j]
    
    gamma_i = 0.5 * (dij + r[i] - r[j])
    gamma_j = dij - gamma_i
    
    return gamma_i, gamma_j

def update_distance_matrix(D, names, i, j):
    """
    Une los nodos i y j en un nuevo nodo k.
    Devuelve:
        - nueva matriz de distancias D'
        - nueva lista de nombres
        - índice del nuevo nodo k
    """
    n = len(D)
    dij = D[i][j]

    # Nombre del nuevo nodo interno
    new_name = f"({names[i]},{names[j]})"

    # Calcular distancias d_km
    new_row = []
    for m in range(n):
        if m != i and m != j:
            dkm = 0.5 * (D[i][m] + D[j][m] - dij)
            new_row.append(dkm)

    # Construir nueva matriz sin i y j
    new_D = []
    new_names = []

    for idx in range(n):
        if idx not in (i, j):
            row = []
            for idx2 in range(n):
                if idx2 not in (i, j):
                    row.append(D[idx][idx2])
            new_D.append(row)
            new_names.append(names[idx])

    # Añadir la columna del nuevo nodo k
    for r in range(len(new_D)):
        new_D[r].append(new_row[r])

    # Añadir la fila del nuevo nodo k
    new_row.append(0.0)
    new_D.append(new_row)

    # Añadir nombre del nuevo nodo
    new_names.append(new_name)

    # Índice del nuevo nodo k
    k_index = len(new_names) - 1

    return new_D, new_names, k_index


def add_to_tree(tree, names, i, j, k_name, gamma_i, gamma_j):
    """
    Añade un nodo interno k al árbol.
    names[i] y names[j] son los hijos.
    """
    tree[k_name] = {
        "left": names[i],
        "right": names[j],
        "length_left": gamma_i,
        "length_right": gamma_j
    }


def build_newick(tree, node):
    """
    Construye el Newick recursivamente desde el nodo dado.
    """
    if node not in tree:
        return node  # es una hoja

    left = tree[node]["left"]
    right = tree[node]["right"]
    L = tree[node]["length_left"]
    R = tree[node]["length_right"]

    left_newick = build_newick(tree, left)
    right_newick = build_newick(tree, right)

    return f"({left_newick}:{L},{right_newick}:{R})"


def neighbor_joining(names, D):
    tree = {}
    current_names = names[:]
    current_D = [row[:] for row in D]

    while len(current_names) > 3:
        # 1. r_i
        r = compute_r_values(current_D)

        # 2. matriz N
        N = compute_N_matrix(current_D, r)

        # 3. par mínimo
        i, j = find_min_N(N)

        # 4. longitudes de rama
        gamma_i, gamma_j = compute_branch_lengths(current_D, r, i, j)

        # 5. nombre del nuevo nodo
        k_name = f"({current_names[i]},{current_names[j]})"

        # 6. registrar en el árbol
        add_to_tree(tree, current_names, i, j, k_name, gamma_i, gamma_j)

        # 7. actualizar matriz de distancias
        current_D, current_names, k_index = update_distance_matrix(
            current_D, current_names, i, j
        )

    # Cuando quedan 3 nodos, conectarlos directamente
    a, b, c = current_names
    # Distancias finales
    dab = current_D[0][1]
    dac = current_D[0][2]
    dbc = current_D[1][2]

    # Fórmulas finales
    La = (dab + dac - dbc) / 2
    Lb = (dab + dbc - dac) / 2
    Lc = (dac + dbc - dab) / 2

    root = f"({a}:{La},{b}:{Lb},{c}:{Lc})"
    return tree, root


def run_nj(filepath):
    names, D = read_phylip_distance_matrix(filepath)
    tree, root = neighbor_joining(names, D)
    newick = build_newick(tree, root) + ";"
    return newick







