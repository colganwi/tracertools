from collections import Counter
import os
import subprocess
import tempfile

import cassiopeia as cas
import networkx as nx
import numpy as np
import pandas as pd
import pycea
import treedata as td

from .config import node_name_generator
from .utils import bfs_names, mask_truncal_edits, newick_to_tree, save_characters_fasta, save_edit_distance


def get_root(tree):
    """Finds the root of a tree"""
    return [node for node in tree.nodes if tree.in_degree(node) == 0][0]


def get_leaves(tree: nx.DiGraph):
    """Finds the leaves of a tree"""
    return [node for node in nx.dfs_postorder_nodes(tree, get_root(tree)) if tree.out_degree(node) == 0]

def split_tree(
    tree: nx.DiGraph,
    threshold: int | None = None,
    root=None,
    node_name_gen: callable = node_name_generator(),
    characters_key: str = "characters",
    root_state="*",
):
    """
    Split a rooted tree into subtrees.

    Default behavior matches the original version: split into one subtree per
    child of the root.

    If threshold is provided, traverse from the root and split at edges where
    n_edits >= threshold.
    """
    if root is None:
        root = get_root(tree)

    subtrees = []

    def _add_artificial_root(subtree, child):
        new_root = next(node_name_gen)

        if characters_key in subtree.nodes[child]:
            n_characters = len(subtree.nodes[child][characters_key])
            subtree.add_node(
                new_root,
                **{characters_key: [root_state] * n_characters},
            )
        else:
            subtree.add_node(new_root)

        subtree.add_edge(new_root, child)
        return subtree

    # Default behavior: split once at children of root
    if threshold is None:
        for child in tree.successors(root):
            nodes = {child} | nx.descendants(tree, child)
            subtree = tree.subgraph(nodes).copy()
            subtree = _add_artificial_root(subtree, child)
            subtrees.append(subtree)

        return subtrees

    # Threshold behavior
    for u, v in tree.edges:
        tree.edges[u, v]["n_edits"] = np.sum(
            np.array(tree.nodes[v][characters_key])
            != np.array(tree.nodes[u][characters_key])
        )

    stack = [root]

    while stack:
        u = stack.pop()

        for v in tree.successors(u):
            n_edits = tree.edges[u, v].get("n_edits", 0)

            if n_edits >= threshold:
                children = list(tree.successors(v))

                c1_descendants = (
                    len(nx.descendants(tree, children[0])) if len(children) == 2 else 0
                )
                c2_descendants = (
                    len(nx.descendants(tree, children[1])) if len(children) == 2 else 0
                )

                if (
                    ((c1_descendants > 20) or (c2_descendants > 20))
                    and ((c1_descendants < 5) or (c2_descendants < 5))
                ):
                    stack.append(v)
                else:
                    nodes = {v} | nx.descendants(tree, v)
                    subtree = tree.subgraph(nodes).copy()
                    subtree = _add_artificial_root(subtree, v)
                    subtrees.append(subtree)
            else:
                stack.append(v)

    subtrees.sort(key=lambda t: t.number_of_nodes(), reverse=True)
    return subtrees


def infer_ancestral_characters(tdata, tree="tree", key="characters", edit_cost=0.6, copy=False):
    """Reconstruct ancestral characters using Sankoff algorithm"""
    tree_key = tree
    if copy:
        tdata = tdata.copy()
    n_characters = max(tdata.obsm[key].max().max() + 1, 10)
    costs = np.ones(shape=(n_characters, n_characters), dtype=float)
    costs[0, :] = edit_cost
    np.fill_diagonal(costs, 0)
    costs = pd.DataFrame(costs, index=range(0, n_characters), columns=range(0, n_characters))
    pycea.tl.ancestral_states(tdata, keys=key, method="sankoff", costs=costs, missing_state=-1, tree=tree_key)
    if copy:
        return tdata


def indices_of_most_common(lst):
    c = Counter(x for x in lst if x not in ("*", "-"))
    if not c:
        return None
    val, count = c.most_common(1)[0]
    return [i for i, x in enumerate(lst) if x == val] if count > 1 else None


def collapse_mutationless_edges(tree, key="characters", mutation_key=None, rescue_edges = False,node_name_gen= node_name_generator(), copy=False):
    """Collapse mutationless edges"""
    if copy:
        tree = tree.copy()
    root = get_root(tree)
    for edge in reversed(list(nx.dfs_edges(tree, root))):
        if mutation_key is not None:
            has_mutation = tree.edges[edge][mutation_key]
        else:
            has_mutation = np.any(tree.nodes[edge[1]][key] != tree.nodes[edge[0]][key])
        if not has_mutation:
            children = list(tree.successors(edge[1]))
            if len(children) > 0:
                for child in children:
                    tree.add_edge(edge[0], child)
                tree.remove_edge(*edge)
                tree.remove_node(edge[1])
    if rescue_edges:
        n_characters = len(tree.nodes[root][key])
        for node in list(tree.nodes):
            children = list(tree.successors(node))
            if len(children) > 2:
                for i in range(n_characters):
                    if tree.nodes[node][key][i] == "*":
                        idx_same_edit = indices_of_most_common(
                            [tree.nodes[child][key][i] for child in children]
                        )
                        if idx_same_edit is not None:
                            # create new internal node
                            new_node = next(node_name_gen)
                            new_characters = tree.nodes[node][key].copy()
                            new_characters[i] = tree.nodes[children[idx_same_edit[0]]][key][i]
                            tree.add_edge(node, new_node)
                            tree.nodes[new_node][key] = new_characters
                            # reassign children
                            for idx in idx_same_edit:
                                child = children[idx]
                                tree.add_edge(new_node, child)
                                tree.remove_edge(node, child)
                            children = list(tree.successors(node))
                            if len(children) <= 2:
                                break
    if copy:
        return tree


def majority_character(characters, min_size=20, min_frac=0.8):
    """Find the majority character in a list of characters"""
    characters = np.array(characters)
    if len(characters) < min_size:
        return -1
    characters = characters[characters != -1]
    if len(characters) == 0:
        return -1
    unique_values, counts = np.unique(characters, return_counts=True)
    for value, count in zip(unique_values, counts, strict=False):
        if count / len(characters) > min_frac:
            return value
    return 0


def same_characters(c1, c2):
    """Check if two arrays are equal, ignoring -1 values"""
    c1 = np.array(c1)
    c2 = np.array(c2)
    mask = (c1 != -1) & (c2 != -1)
    if not mask.any():
        return True
    return np.array_equal(c1[mask], c2[mask])


def identify_mutations(tdata, tree="tree", key="characters", key_added="has_mutation", min_frac=0.75, copy=False):
    """Mark edges with a mutation"""
    if copy:
        tdata = tdata.copy()
    method = lambda x: majority_character(x, min_frac=min_frac)
    pycea.tl.ancestral_states(tdata, keys=key, method=method, keys_added="majority_characters", tree=tree)
    for edge in tdata.obst[tree].edges:
        has_mutation = not (
            same_characters(
                tdata.obst[tree].nodes[edge[0]]["majority_characters"],
                tdata.obst[tree].nodes[edge[1]]["majority_characters"],
            )
            and same_characters(tdata.obst[tree].nodes[edge[0]][key], tdata.obst[tree].nodes[edge[1]][key])
        )
        tdata.obst[tree].edges[edge][key_added] = has_mutation
    for node in tdata.obst[tree].nodes:
        del tdata.obst[tree].nodes[node]["majority_characters"]
    if copy:
        return tdata

def estimate_leaf_fitness(tdata, tree="tree", depth_key="depth", key_added="fitness", copy=False):
    """Estimate leaf fitness based on tree and branch lengths"""
    tree_key = tree
    if copy:
        tdata = tdata.copy()
    nx_tree = tdata.obst[tree_key].copy()
    for node in nx_tree:
        nx_tree.nodes[node]["_depth"] = nx_tree.nodes[node][depth_key]
    cas_tree = cas.data.CassiopeiaTree(tree=nx_tree)
    for edge in cas_tree.depth_first_traverse_edges():
        t1 = cas_tree.get_attribute(edge[0], "_depth")
        t2 = cas_tree.get_attribute(edge[1], "_depth")
        cas_tree.set_branch_length(edge[0], edge[1], abs(t1 - t2))
    fitness_estimator = cas.tools.fitness_estimator.LBIJungle()
    fitness_estimator.estimate_fitness(cas_tree)
    fitnesses = np.array([cas_tree.get_attribute(cell, "fitness") for cell in cas_tree.leaves])
    fitnesses = pd.Series(fitnesses, index=cas_tree.leaves)
    tdata.obs[key_added] = tdata.obs_names.map(fitnesses)
    if copy:
        return tdata


def estimate_branch_lengths(
    tree: nx.DiGraph,
    key: str = "characters",
    method: str = "convexml",
    minimum_branch_length=0.001,
    key_added: str = "time",
    solver: str = "CLARABEL",
    verbose: bool = False,
    missing_state = "-",
    unedited_state = "*",
    mutation_rates = None,
    pseudo_count = 0.5,
    total_time = 1.0,
    lamlpro_cmd: str = "lamlpro",
):
    """Estimate branch lengths for a tree.

    Parameters
    ----------
    tree : the tree to estimate branch lengths for
    key : the node attribute holding the character states
    method : 'convexml' (IIDExponentialMLE, default) or 'laml' (LAML-Pro optimize mode)
    minimum_branch_length : minimum branch length (convexml only)
    key_added : the node attribute to store the estimated times under
    solver : the solver to use for IIDExponentialMLE (convexml only)
    verbose : whether to print verbose output
    missing_state, unedited_state : encoding of missing and unedited states
    mutation_rates : relative mutation rates (convexml only)
    pseudo_count : pseudo mutations per edge (convexml only)
    total_time : root-to-leaf time the tree is scaled to (None to leave unscaled)
    lamlpro_cmd : path to / name of the ``lamlpro`` binary (laml only)
    """
    leaves = get_leaves(tree)
    if method == "convexml":
        from convexml._iid_exponential_mle import IIDExponentialMLE
        chr_to_int = {str(i):i for i in range(9)}
        chr_to_int.update({missing_state:-1, unedited_state:0, "!":9})
        node_states = {}
        for node in tree.nodes:
            char_states = tree.nodes[node][key]
            node_states[node] = [chr_to_int[char] for char in char_states]
        characters = pd.DataFrame(node_states).T.loc[leaves]
        cas_tree = cas.data.CassiopeiaTree(character_matrix=characters, tree=tree)
        cas_tree.set_all_character_states(node_states)
        if mutation_rates is not None:
            mutation_rates = mutation_rates * int(characters.shape[1]/len(mutation_rates))
        IIDExponentialMLE(
            verbose=verbose, solver=solver, minimum_branch_length=minimum_branch_length, relative_mutation_rates=mutation_rates, pseudo_mutations_per_edge = pseudo_count
        ).estimate_branch_lengths(cas_tree)
        node_times = nx.get_node_attributes(cas_tree.get_tree_topology(), "time")
        if total_time is not None:
            for node in node_times:
                node_times[node] = node_times[node] * total_time
    elif method == "laml":
        from .solver import laml
        characters = pd.DataFrame(
            {leaf: list(tree.nodes[leaf][key]) for leaf in leaves}
        ).T
        laml_tree = laml(
            characters,
            tree,
            mode="optimize",
            ultrametric=True,
            missing_state=missing_state,
            unedited_state=unedited_state,
            lamlpro_cmd=lamlpro_cmd,
            min_branch_length=minimum_branch_length,
        )
        node_times = nx.get_node_attributes(laml_tree, "time")
        if total_time is not None:
            height = max(node_times.values())
            if height > 0:
                scale = total_time / height
                for node in node_times:
                    node_times[node] = node_times[node] * scale
        # A root unifurcation (root -> single child) is collapsed inside laml, so the
        # original root is absent from node_times. It marks time 0 (the founder), and
        # its child's time is the LAML-estimated time before the first division.
        root = get_root(tree)
        if root not in node_times:
            node_times[root] = 0.0
    else:
        raise ValueError(f"Unknown method: {method}. Use 'convexml' or 'laml'.")
    nx.set_node_attributes(tree, node_times, key_added)


def collapse_unifurcations(G: nx.DiGraph, allow_root: bool = True) -> nx.DiGraph:
    """Collapse unifurcations in a tree.

    Contract every node with out_degree == 1 by bypassing it.
    If allow_root=True, also removes the root when it has a single child.
    Edge attributes from the bypassed edge (u->v) are preserved.
    """
    changed = True
    while changed:
        changed = False
        to_contract = [
            n for n in list(G.nodes())
            if G.out_degree(n) == 1 and (G.in_degree(n) == 1 or (allow_root and G.in_degree(n) == 0))
        ]
        if not to_contract:
            break
        for n in to_contract:
            # Node may have been removed in this round
            if n not in G:
                continue
            succ = next(G.successors(n))  # unique child
            preds = list(G.predecessors(n))  # 0 or 1 in a tree
            edge_attrs = G.get_edge_data(n, succ) or {}
            if preds:  # normal degree-2 case
                p = preds[0]
                if p != succ:
                    G.add_edge(p, succ, **edge_attrs)
            # Remove n (and its incident edges)
            if n in G:
                G.remove_node(n)
            changed = True
    return G

def replace_subtrees(tree: nx.DiGraph, replacements: list[nx.DiGraph], error_on_missing = False) -> nx.DiGraph:
    """
    Replace subtrees of `tree` with the given replacement trees.

    Each replacement DiGraph must be a rooted tree whose root node
    identifier exists in `tree`. For each replacement R:

      * Let r be the root of R (in_degree == 0 in R).
      * In `tree`, remove all descendants of r (but keep r itself).
      * Graft in R under r: update r's node attributes from R, add
        new nodes/edges from R.

    The incoming edge into r from its parent in `tree` is preserved.
    Returns a new DiGraph; does not modify `tree` in place.
    """
    tree = tree.copy()

    for replacement in replacements:
        # --- find root of replacement tree ---
        root = get_root(replacement)

        if root not in tree:
            if error_on_missing:
                raise KeyError(f"Replacement root {root!r} not found in base tree")
            else:
                continue

        # --- remove old subtree below r (but keep r itself) ---
        descendants = nx.descendants(tree, root)
        tree.remove_nodes_from(descendants)

        # At this point r has no outgoing edges.

        # --- add / merge nodes from replacement ---
        for n, data in replacement.nodes(data=True):
            if n == root:
                # Update attributes on existing root
                tree.nodes[root].update(data)
            else:
                if n in tree:
                    raise ValueError(
                        f"Node {n!r} from replacement subtree already exists in base tree"
                    )
                tree.add_node(n, **data)
        # --- add edges from replacement ---
        for u, v, data in replacement.edges(data=True):
            # In R, root r has no incoming edges, so we only add its outgoing edges
            tree.add_edge(u, v, **data)
    return tree

def ancestral_characters(
    tree: nx.DiGraph,
    characters: pd.DataFrame,
    missing_state: str = "-",
    unedited_state: str = "*",
    reedit_penalty: float = 2.5,
    key_added: str = "characters",
) -> None:
    """Reconstruct ancestral characters using weighted Sankoff algorithm.

    Parameters
    ----------
    tree
        Directed tree (root has in-degree 0).
    characters
        DataFrame of leaf characters; index should match leaf node ids,
        columns are character positions.
    missing_state
        Symbol used for missing data.
    unedited_state
        Special baseline state.
    reedit_penalty
        Cost of re-editing away from the unedited state.
    key_added
        Node attribute key where reconstructed characters are stored.
    """
    alphabet, val_to_idx, cost_matrix = _build_cost_matrix(
        characters, missing_state, unedited_state, reedit_penalty
    )

    n_chars = characters.shape[1]
    n_states = len(alphabet)

    # Initialize all nodes with None
    node_attrs = {node: [None] * n_chars for node in tree.nodes}
    for node, row in characters.iterrows():
        if node in node_attrs:
            node_attrs[node] = row.to_list()
    nx.set_node_attributes(tree, node_attrs, key_added)

    _reconstruct_sankoff_weighted(
        tree=tree,
        key=key_added,
        alphabet=alphabet,
        val_to_idx=val_to_idx,
        cost_matrix=cost_matrix,
        missing=missing_state,
        default=unedited_state,
    )


def _build_cost_matrix(characters, missing_state, unedited_state, reedit_penalty):
    vals = pd.unique(characters.to_numpy().ravel())
    vals = vals[(vals != missing_state) & (vals != unedited_state)]
    alphabet = [unedited_state, *vals.tolist()]
    n_states = len(alphabet)

    M = np.full((n_states, n_states), reedit_penalty, float)
    M[0, :] = 1.0
    M[:, 0] = reedit_penalty
    np.fill_diagonal(M, 0.0)

    val_to_idx = {v: i for i, v in enumerate(alphabet)}
    return alphabet, val_to_idx, M

def _reconstruct_sankoff_weighted(
    tree: nx.DiGraph,
    key: str,
    alphabet: list[str],
    val_to_idx: dict[str, int],
    cost_matrix: np.ndarray,
    missing: str,
    default: str,
):
    n_states = len(alphabet)
    n_chars = len(next(iter(tree.nodes(data=True)))[1][key])

    default_idx = val_to_idx.get(default, None)
    children = {n: list(tree.successors(n)) for n in tree.nodes}
    root = next(n for n, d in tree.in_degree() if d == 0)

    postorder = list(nx.dfs_postorder_nodes(tree, root))

    # DP matrices
    scores = {}       # node -> (n_states, n_chars)
    sizes = {}        # node -> subtree size
    pointers = {}     # node -> (n_states, n_children, n_chars)

    # Upward pass -------------------------------------------------------------

    for node in postorder:
        ch = children[node]

        # LEAF
        if not ch:
            leaf_vals = np.array(tree.nodes[node][key])
            leaf_missing = (leaf_vals == missing)

            # Inflate to (n_states × n_chars)
            S = np.full((n_states, n_chars), np.inf)

            for c, v in enumerate(leaf_vals):
                if leaf_missing[c]:
                    S[:, c] = 0
                else:
                    S[val_to_idx[v], c] = 0

            scores[node] = S
            sizes[node] = 1
            continue

        # INTERNAL NODE
        num_children = len(ch)
        node_scores = np.zeros((n_states, n_chars))
        ptrs = np.empty((n_states, num_children, n_chars), dtype=np.int32)
        total_size = 0

        # Vectorized DP for each child
        for i, child in enumerate(ch):
            c_scores = scores[child]          # (n_states, n_chars)
            c_size = sizes[child]
            total_size += c_size
            weight = np.log2(1 + c_size)

            # Weighted cost: parent p, child state s
            # cost_matrix[p,s] shape => (n_states, n_states)
            # Expand to (n_states, n_states, n_chars)
            weighted = c_scores[None, :, :] + weight * cost_matrix[:, :, None]

            # Choose s that minimizes cost for each (p, char)
            best_idx = np.argmin(weighted, axis=1)              # (n_states, n_chars)
            best_cost = weighted[np.arange(n_states)[:,None], best_idx, np.arange(n_chars)[None,:]]

            node_scores += best_cost
            ptrs[:, i, :] = best_idx

        scores[node] = node_scores
        sizes[node] = total_size
        pointers[node] = ptrs

    # Downward pass -----------------------------------------------------------

    root_scores = scores[root]

    if default_idx is not None:
        root_assign = np.full(n_chars, default_idx, dtype=int)
    else:
        root_assign = np.argmin(root_scores, axis=0).astype(int)

    # Assign root
    root_states = [alphabet[i] for i in root_assign]
    tree.nodes[root][key] = root_states

    stack = [(root, root_assign)]

    while stack:
        node, parent_states = stack.pop()
        if node not in pointers:
            continue

        ptrs = pointers[node]   # (n_states, n_children, n_chars)

        for child_idx, child in enumerate(children[node]):
            child_assign = ptrs[parent_states, child_idx, np.arange(n_chars)]
            tree.nodes[child][key] = [alphabet[i] for i in child_assign]
            stack.append((child, child_assign))

def count_branch_edits(tree, key="characters"):
    for u, v in tree.edges():
        n = sum(a != b for a, b in zip(tree.nodes[u][key], tree.nodes[v][key]))
        tree.edges[u, v]["n_edits"] = n

