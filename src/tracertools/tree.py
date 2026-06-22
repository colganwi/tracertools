import os
import subprocess
import tempfile
from collections import Counter

import cassiopeia as cas
import networkx as nx
import numpy as np
import pandas as pd
import pycea

from .config import node_name_generator
from .utils import bfs_names, mask_truncal_edits, newick_to_tree, save_characters_fasta, save_edit_distance


def get_root(tree):
    """Finds the root of a tree"""
    return [node for node in tree.nodes if tree.in_degree(node) == 0][0]


def get_leaves(tree: nx.DiGraph):
    """Finds the leaves of a tree"""
    return [node for node in nx.dfs_postorder_nodes(tree, get_root(tree)) if tree.out_degree(node) == 0]


def split_tree(tree: nx.DiGraph, node_name_gen: callable = node_name_generator()):
    """
    Split a rooted tree into subtrees, one per child of the root.

    Parameters
    ----------
    tree : nx.DiGraph
        A directed tree with a single root (in-degree 0).
    node_name_gen : callable
        A generator function for naming nodes in the subtrees.

    Returns
    -------
    list[nx.DiGraph]
        List of DiGraphs, each being the subtree rooted at one child of the root.
    """
    # find root (node with in-degree 0)
    root = get_root(tree)

    subtrees = []
    for child in tree.successors(root):
        nodes = {child} | nx.descendants(tree, child)
        subtree = tree.subgraph(nodes).copy()
        subtree.add_edge(next(node_name_gen), child)
        subtrees.append(subtree)
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


def collapse_mutationless_edges(
    tree, key="characters", mutation_key=None, rescue_edges=False, node_name_gen=node_name_generator(), copy=False
):
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
                        idx_same_edit = indices_of_most_common([tree.nodes[child][key][i] for child in children])
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


def reconstruct_fasttree(
    character_matrix: pd.DataFrame,
    mutation_rate: float = 0.1,
    number_of_states: int = 10,
    add_root: bool = False,
    fasttree_cmd: str = "FastTree",
    logfile: str | None = None,
) -> nx.DiGraph:
    """Reconstruct a lineage tree with FastTree.

    Parameters
    ----------
    character_matrix : pd.DataFrame
        Rows = samples (index), cols = sites, values in 0..number_of_states-1 or -1.
    mutation_rate : float
    number_of_states : int (<= 20)
    add_root : bool
        If True, midpoint-root the ete3 tree (using set_outgroup(get_midpoint_outgroup())).
    fasttree_cmd : str
    logfile : Optional[str]

    Returns
    -------
    nx.DiGraph
    """
    # validation
    if not (0 <= mutation_rate <= 1):
        raise ValueError("mutation_rate must be between 0 and 1.")
    if number_of_states > 20:
        raise ValueError("number_of_states must be <= 20.")
    if character_matrix.isnull().values.any():
        raise ValueError("character_matrix contains NaNs; use -1 for missing.")

    vals = character_matrix.values
    if not np.isin(vals[vals != -1], np.arange(number_of_states)).all():
        raise ValueError("character_matrix has values outside 0..number_of_states-1 (excluding -1).")

    detected = np.unique(vals[vals != -1])
    if len(detected) != number_of_states:
        raise ValueError(f"Detected {len(detected)} states, but number_of_states={number_of_states}.")

    with tempfile.TemporaryDirectory() as td:
        fasta_path = os.path.join(td, "character_matrix.fasta")
        edit_base = os.path.join(td, "edit_distance")
        tree_path = os.path.join(td, "tree.nwck")

        save_edit_distance(number_of_states, edit_base)
        save_characters_fasta(character_matrix, fasta_path)

        cmd = (
            f"bash -c '. ~/.bashrc && "
            f"{fasttree_cmd} -nosupport -noml -nome "
            f"-matrix {edit_base} "
            f"{fasta_path} > {tree_path}'"
        )

        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if logfile:
            with open(logfile, "ab") as lf:
                lf.write(stdout or b"")
                lf.write(stderr or b"")
        if process.returncode != 0:
            raise RuntimeError(f"FastTree failed (code {process.returncode}). Stderr:\n{stderr.decode('utf-8')}")
        # Read and (optionally) root with ete3
        with open(tree_path) as f:
            newick_str = f.read().strip()
            tree = newick_to_tree(newick_str, midpoint_root=True if add_root else False)

    return tree


def reconstruct_cassiopeia(characters, solver="upgma"):
    """Reconstruct a tree using cassiopeia"""
    solvers = {
        "nj": cas.solver.NeighborJoiningSolver(
            cas.solver.dissimilarity.weighted_hamming_distance, add_root=True, fast=True
        ),
        "upgma": cas.solver.UPGMASolver(cas.solver.dissimilarity.weighted_hamming_distance, fast=True),
        "greedy": cas.solver.VanillaGreedySolver(),
    }
    cas_tree = cas.data.CassiopeiaTree(character_matrix=characters)
    solvers[solver].solve(cas_tree)
    tree = cas_tree.get_tree_topology()
    tree = nx.relabel_nodes(tree, bfs_names(tree))
    return tree


def reconstruct_tree(
    tdata, solver="fasttree", key="characters", tree_added="tree", mask_truncal=True, logfile=None, copy=False
):
    """Reconstruct a tree from characters optionally estimating branch lengths"""
    if copy:
        tdata = tdata.copy()
    characters = tdata.obsm[key]
    if mask_truncal:
        characters = mask_truncal_edits(characters).copy()
    if solver == "fasttree":
        tree = reconstruct_fasttree(characters, add_root=True, logfile=logfile)
    else:
        tree = reconstruct_cassiopeia(characters, solver=solver)
    tdata.obst[tree_added] = tree
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
    minimum_branch_length=0.001,
    key_added: str = "time",
    solver: str = "CLARABEL",
    verbose: bool = False,
    missing_state="-",
    unedited_state="*",
    mutation_rates=None,
    pseudo_count=0.5,
    total_time=1.0,
):
    """Estimate branch lengths using IIDExponentialMLE

    Parameters
    ----------
    tree : the tree to estimate branch lengths for
    key : the key for the character matrix in tdata.obsm
    key_added : the key for the added branch lengths in tdata.obst
    solver : the solver to use for IIDExponentialMLE
    verbose : whether to print verbose output
    copy : whether to copy the TreeData object

    Returns
    -------
    tdata : the modified TreeData object if copy is True
    """
    from convexml._iid_exponential_mle import IIDExponentialMLE

    leaves = get_leaves(tree)
    chr_to_int = {str(i): i for i in range(9)}
    chr_to_int.update({missing_state: -1, unedited_state: 0, "!": 9})
    node_states = {}
    for node in tree.nodes:
        char_states = tree.nodes[node][key]
        node_states[node] = [chr_to_int[char] for char in char_states]
    characters = pd.DataFrame(node_states).T.loc[leaves]
    cas_tree = cas.data.CassiopeiaTree(character_matrix=characters, tree=tree)
    cas_tree.set_all_character_states(node_states)
    if mutation_rates is not None:
        mutation_rates = mutation_rates * int(characters.shape[1] / len(mutation_rates))
    IIDExponentialMLE(
        verbose=verbose,
        solver=solver,
        minimum_branch_length=minimum_branch_length,
        relative_mutation_rates=mutation_rates,
        pseudo_mutations_per_edge=pseudo_count,
    ).estimate_branch_lengths(cas_tree)
    node_times = nx.get_node_attributes(cas_tree.get_tree_topology(), "time")
    if total_time is not None:
        for node in node_times:
            node_times[node] = node_times[node] * total_time
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
            n
            for n in list(G.nodes())
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


def replace_subtrees(tree: nx.DiGraph, replacements: list[nx.DiGraph], error_on_missing=False) -> nx.DiGraph:
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
                    raise ValueError(f"Node {n!r} from replacement subtree already exists in base tree")
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
    alphabet, val_to_idx, cost_matrix = _build_cost_matrix(characters, missing_state, unedited_state, reedit_penalty)

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
    scores = {}  # node -> (n_states, n_chars)
    sizes = {}  # node -> subtree size
    pointers = {}  # node -> (n_states, n_children, n_chars)

    # Upward pass -------------------------------------------------------------

    for node in postorder:
        ch = children[node]

        # LEAF
        if not ch:
            leaf_vals = np.array(tree.nodes[node][key])
            leaf_missing = leaf_vals == missing

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
            c_scores = scores[child]  # (n_states, n_chars)
            c_size = sizes[child]
            total_size += c_size
            weight = np.log2(1 + c_size)

            # Weighted cost: parent p, child state s
            # cost_matrix[p,s] shape => (n_states, n_states)
            # Expand to (n_states, n_states, n_chars)
            weighted = c_scores[None, :, :] + weight * cost_matrix[:, :, None]

            # Choose s that minimizes cost for each (p, char)
            best_idx = np.argmin(weighted, axis=1)  # (n_states, n_chars)
            best_cost = weighted[np.arange(n_states)[:, None], best_idx, np.arange(n_chars)[None, :]]

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

        ptrs = pointers[node]  # (n_states, n_children, n_chars)

        for child_idx, child in enumerate(children[node]):
            child_assign = ptrs[parent_states, child_idx, np.arange(n_chars)]
            tree.nodes[child][key] = [alphabet[i] for i in child_assign]
            stack.append((child, child_assign))


def count_branch_edits(tree, key="characters"):
    for u, v in tree.edges():
        n = sum(a != b for a, b in zip(tree.nodes[u][key], tree.nodes[v][key]))
        tree.edges[u, v]["n_edits"] = n
