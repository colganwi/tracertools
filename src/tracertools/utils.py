from typing import Any, Dict, Optional

import ete3
import networkx as nx
import numpy as np
import pandas as pd
from numba import njit
from pycea.utils import get_leaves, get_root

NUM_TO_AA = {
    0:'A',1:'R',2:'N',3:'D',4:'C',5:'Q',6:'E',7:'G',8:'H',
    9:'I',10:'L',11:'K',12:'M',13:'F',14:'P',15:'S',16:'T',17:'W',18:'Y',19:'V',-1:'-'
}
AAS = [NUM_TO_AA[i] for i in range(20)]

def save_edit_distance(number_of_states: int, basepath: str) -> None:
    """Generate a random edit distance matrix with eigen decomposition.

    Parameters
    ----------
    number_of_states : int
        Number of states to initialize with specific values in the distance matrix.
    basepath : str
        Base path (without extension) to save the output files:
          - `<basepath>.distances`   : distance matrix
          - `<basepath>.inverses`    : inverse eigenvectors
          - `<basepath>.eigenvalues` : eigenvalues

    Notes
    -----
    - The function does not return anything; results are written to files.
    - Requires `AAS`, a list of 20 amino acid symbols, to be defined globally.
    """
    # Step 1: Generate symmetric random matrix
    distances = np.random.uniform(-0.5, 0.5, (20, 20))
    distances = distances @ distances.T

    # Step 2: Apply constraints to the top-left [number_of_states x number_of_states] block
    used = np.ones((number_of_states, number_of_states)) * 2
    used[:, 0] = 1
    used[0, :] = 1
    np.fill_diagonal(used, 0)
    distances[0:number_of_states, 0:number_of_states] = used

    # Step 3: Eigen decomposition
    eigenval, eigenvectors = np.linalg.eig(distances)
    eigeniv = np.real(np.linalg.inv(eigenvectors))
    eigenval = np.real(eigenval)

    # Step 4: Save results to files
    with open(basepath + ".distances", "w") as f:
        f.write("\t".join(AAS) + "\n")
        for i in range(20):
            f.write(AAS[i])
            for j in range(20):
                f.write(f"\t{distances[i, j]}")
            f.write("\n")

    with open(basepath + ".inverses", "w") as f:
        f.write("\t".join(AAS) + "\n")
        for i in range(20):
            f.write(AAS[i])
            for j in range(20):
                f.write(f"\t{eigeniv[i, j]}")
            f.write("\n")

    with open(basepath + ".eigenvalues", "w") as f:
        for i in range(20):
            f.write(f"{eigenval[i]}\n")
        f.write("\n")


def save_characters_fasta(character_matrix: pd.DataFrame, path: str) -> None:
    """Save the character matrix as a FASTA file.

    Parameters
    ----------
    character_matrix : pd.DataFrame
        The character matrix to save.
    path : str
        The path to the output FASTA file.
    """
    with open(path, "w") as f:
        for i in range(character_matrix.shape[0]):
            aa_seq = "".join(character_matrix.iloc[i].map(NUM_TO_AA))
            f.write(f">{character_matrix.index[i]}\n{aa_seq}\n")


def tree_to_newick(
    tree: nx.DiGraph,
    record_branch_lengths: bool = False,
    record_node_names: bool = False,
) -> str:
    """Converts a networkx graph to a newick string.

    Parameters
    ----------
        tree: A networkx tree
        record_branch_lengths: Whether to record branch lengths on the tree in
            the newick string
        record_node_names: Whether to record internal node names on the tree in
            the newick string

    Returns
    -------
        A newick string representing the topology of the tree
    """
    def _to_newick_str(g, node):
        is_leaf = g.out_degree(node) == 0
        weight_string = ""
        if record_branch_lengths and g.in_degree(node) > 0:
            parent = list(g.predecessors(node))[0]
            weight_string = ":" + str(g[parent][node]["length"])
        _name = str(node)
        name_string = ""
        if record_node_names:
            name_string = f"{_name}"
        return (
            "%s" % (_name,) + weight_string
            if is_leaf
            else (
                "("
                + ",".join(
                    _to_newick_str(g, child) for child in g.successors(node)
                )
                + ")"
                + name_string
                + weight_string
            )
        )
    root = [node for node in tree if tree.in_degree(node) == 0][0]
    return _to_newick_str(tree, root) + ";"


def _safe_label(name: str) -> str:
    """Quote labels that contain Newick-reserved characters."""
    s = str(name)
    if any(ch in s for ch in "():,; \t'"):
        # Newick escapes single quotes by doubling them
        s = "'" + s.replace("'", "''") + "'"
    return s


def newick_to_tree(
    newick: str,
    length_attr: Optional[str] = "length",
    midpoint_root: bool = False,
) -> nx.DiGraph:
    """Parse Newick via ete3 and return a directed NetworkX DiGraph).

    Parameters
    ----------
    newick : str
        Newick string (with trailing ';').
    length_attr : str, optional
        Edge attribute name to store branch lengths (child.dist). If None, lengths are dropped.
    midpoint_root : bool, default False
        If True, root the ete3 tree at its midpoint before conversion.

    Returns
    -------
    nx.DiGraph
    """
    T = ete3.Tree(newick, format=1)  # robust default for standard Newick
    if midpoint_root:
        T.set_outgroup(T.get_midpoint_outgroup())
    G = nx.DiGraph()
    counter = 0

    def node_name(n: ete3.Tree) -> str:
        nonlocal counter
        if n.name and n.name.strip():
            return n.name
        counter += 1
        return f"node{counter}"

    # assign stable names
    for n in T.traverse("preorder"):
        n.add_feature("_nx_name", node_name(n))

    # add edges parent->child, with optional branch lengths from child.dist
    for n in T.traverse("preorder"):
        parent = n._nx_name
        G.add_node(parent)
        for c in n.children:
            child = c._nx_name
            if length_attr is not None:
                G.add_edge(parent, child, **{length_attr: float(getattr(c, "dist", 0.0) or 0.0)})
            else:
                G.add_edge(parent, child)

    return G

@njit
def hamming_distance(arr1, arr2):
    """Compute Hamming distance between two cells character arrays"""
    valid_mask = (arr1 != -1) & (arr2 != -1)
    hamming_distance = 0
    for x, y in zip(arr1[valid_mask], arr2[valid_mask]):
        if x == y:
            pass
        elif x == 0 or y == 0:
            hamming_distance += 1
        else:
            hamming_distance += 2
    num_valid_comparisons = np.sum(valid_mask)
    if num_valid_comparisons == 0:
        return 0
    normalized_distance = hamming_distance / num_valid_comparisons
    return normalized_distance


def bfs_names(tree):
    """Assign names to nodes in BFS order"""
    bfs_nodes = list(nx.bfs_tree(tree, source=get_root(tree)))
    leaves = get_leaves(tree)
    return {node: f"node{i}" for i, node in enumerate(bfs_nodes) if node not in leaves}


def get_edit_frac(characters):
    """Compute the fraction of edited characters."""
    characters = np.array(characters)
    n_detected = np.sum(characters > -1,axis = 1)
    n_edited = np.sum(characters > 0,axis = 1)
    return n_edited/n_detected


def get_detection_rate(characters):
    """Compute the fraction of detected characters."""
    return np.sum(characters > 0)/len(characters)


def mask_truncal_edits(characters):
    """Mask characters with edit shared by all cells"""
    masked_characters = {}
    for column in characters.columns:
        filtered_values = characters[characters[column] != -1][column]
        value_counts = filtered_values.value_counts()
        if len(value_counts[(value_counts.index != 0) & (value_counts.index != -1)]) > 0:
            most_common_value = value_counts[(value_counts.index != 0) & (value_counts.index != -1)].idxmax()
            fraction = (filtered_values == most_common_value).sum() / len(filtered_values)
            if fraction < .95:
                masked_characters[column] = characters[column]
    return pd.DataFrame(masked_characters)