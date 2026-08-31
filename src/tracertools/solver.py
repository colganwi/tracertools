import os
import subprocess
import tempfile

import ete3
import networkx as nx
import numpy as np
import pandas as pd

from .config import node_name_generator
from .root import centroid, character_outgroup
from .tree import collapse_unifurcations
from .utils import (
    get_leaves,
    get_root,
    newick_to_tree,
    save_characters_csv,
    save_characters_fasta,
    save_edit_distance,
    tree_to_newick,
)

Mutation = tuple[str, str]  # (column, value)


def _connect_leaves(tree):
    new = tree.copy()
    for node, data in tree.nodes(data=True):
        if tree.out_degree(node) == 0:
            for child in data["descendants"]:
                if data.get("doublet", False):
                    new.add_node(child, depth=data.get("depth", 0) + 1, parent=node, doublet=True)
                new.add_node(child, depth=data.get("depth", 0) + 1, parent=node)
                new.add_edge(node, child)
    return new


def n_mutation_greedy(
    characters: pd.DataFrame,
    jaccard_threshold: float = 0.9,
    overlap_threshold: float = 0.9,
    doublet_threshold: float = 0.1,
    min_clade_frac: float = 0.01,
    max_depth: int = 3,
    min_mutations: int = 3,
    node_name_gen: callable = node_name_generator(),
    connect_leaves: bool = True,
) -> nx.DiGraph:
    """
    Greedily construct a phylogenetic tree using splits defined by >= min_mutations.

    Parameters
    ----------
    characters : pd.DataFrame
        Shape (n_cells, n_characters). Values:
        "*" = no mutation, "-" = missing, anything else = specific mutation state.
        Index labels are treated as sequence IDs.
    jaccard_threshold : float, optional
        Minimum Jaccard similarity required between mutation presence sets
        (vs the seed mutation) to include a mutation in a split.
    overlap_threshold : float, optional
        Minimum fraction of cells in a clade that must have a mutation
        (counting '-' as a match) to consider it truncal in that clade.
    doublet_threshold : float, optional
        Fraction of total cells below which a residual clade is marked as a doublet.
    min_clade_frac : float, optional
        Minimum fraction of total cells required for a clade.
    min_mutations : int, optional
        Minimum number of mutations required for a split.
    node_name_gen : Callable, optional
        Generator function that yields unique node IDs.
    connect_leaves : bool, optional
        If True, add edges from internal nodes to their leaf descendants.

    Returns
    -------
    tree : nx.DiGraph
        Directed graph of clades. Nodes are named by node_name_gen and have attributes:
        - "descendants": list of cells IDs in that clade.
        - "doublet": bool, True if this clade contains doublets.
    doublets : List[np.ndarray]
        List of arrays of cell IDs marked as doublets.
    """
    # Setup
    cell_ids = characters.index.to_numpy()
    col_ids = characters.columns.to_numpy()
    characters = characters.to_numpy(dtype=str)  # shape (n_cell, n_char)
    n_cell, n_char = characters.shape
    min_clade_size = int(n_cell * min_clade_frac)
    col_to_idx: dict[str, int] = {c: i for i, c in enumerate(col_ids)}
    tree = nx.DiGraph()

    # All sets in helpers below are sets of row indices (ints)
    def _mutation_set(
        parent_idx: np.ndarray,
        mutation: Mutation,
    ) -> set[int]:
        """Set of row indices with this (col, value) (exact match, no '-' special-casing)."""
        col, val = mutation
        j = col_to_idx[col]
        col_vals = characters[parent_idx, j]
        mask = col_vals == val
        return set(parent_idx[mask])

    def _mutation_set_with_missing(
        parent_idx: np.ndarray,
        mutation: Mutation,
    ) -> set[int]:
        """Set of row indices that match this (col, value), where '-' counts as a match."""
        col, val = mutation
        j = col_to_idx[col]
        col_vals = characters[parent_idx, j]
        mask = (col_vals == val) | (col_vals == "-")
        return set(parent_idx[mask])

    def _jaccard(a: set[int], b: set[int]) -> float:
        """Jaccard similarity between two sets."""
        if not a and not b:
            return 0.0
        union_size = len(a | b)
        if union_size == 0:
            return 0.0
        return len(a & b) / union_size

    def _enumerate_mutations(
        parent_idx: np.ndarray,
        truncal_muts: set[Mutation],
    ) -> list[Mutation]:
        """All (column, value) pairs present in this clade, excluding '*'/'-'.

        and anything already used in this subtree.
        """
        muts: list[Mutation] = []
        # Work column-wise using NumPy
        for j in range(n_char):
            col_label = col_ids[j]
            col_vals = characters[parent_idx, j]
            unique_vals = np.unique(col_vals)
            for v in unique_vals:
                if v in {"*", "-"}:
                    continue
                m: Mutation = (col_label, str(v))
                if m in truncal_muts:
                    continue
                muts.append(m)
        return muts

    def _match_set_with_missing(
        parent_idx: np.ndarray,
        mutations: list[Mutation],
        max_mismatch_frac: float = 0.1,
    ) -> set[int]:
        """
        Row indices that match all but at most ONE mutation, where '-' counts as a match.

        For each (col, val), a row is considered matching that mutation if
        df[col] == val or df[col] == '-'. Rows are kept if they have
        mismatches for at most one of the mutations.
        """
        if not mutations:
            return set(parent_idx)

        max_mismatches = int(max_mismatch_frac * len(mutations))
        n = len(parent_idx)
        mismatches = np.zeros(n, dtype=np.int8)

        for col, val in mutations:
            j = col_to_idx[col]
            col_vals = characters[parent_idx, j]
            matches = (col_vals == val) | (col_vals == "-")
            mismatches += ~matches

        mask = mismatches <= max_mismatches
        return parent_idx[mask]

    def _high_prevalence_mutations(
        parent_idx: np.ndarray,
        threshold: float = 0.95,
    ) -> set[Mutation]:
        """Get high-prevalence mutations in this clade.

        Return (column, value) mutations where the *most common* non-*,- value
        in that column is present in > threshold of rows in parent_idx,
        counting '-' as a match.
        """
        result: set[Mutation] = set()
        n = len(parent_idx)
        if n == 0:
            return result

        cutoff = threshold * n

        for j in range(n_char):
            col_label = col_ids[j]
            col_vals = characters[parent_idx, j]

            # Exclude '*' and '-' before finding the most common value
            valid_mask = (col_vals != "*") & (col_vals != "-")
            if not valid_mask.any():
                continue

            valid_vals = col_vals[valid_mask]
            unique_vals, counts = np.unique(valid_vals, return_counts=True)

            # Most common mutation value in this column
            best_idx = int(counts.argmax())
            best_val = unique_vals[best_idx]
            best_count = int(counts[best_idx])

            # '-' counts as a match
            dash_count = int((col_vals == "-").sum())

            if (best_count + dash_count > cutoff) & (best_count > cutoff * 0.95) & (best_val != "9"):
                result.add((col_label, str(best_val)))

        return result

    def _find_best_split(
        parent_idx: np.ndarray,
        truncal_muts: set[Mutation],
        min_mutations: int = min_mutations,
    ) -> tuple[list[Mutation], set[int], set[Mutation]]:
        """
        Find one greedy split in this clade:

        - Choose a seed mutation (most prevalent).
        - Add n mutations with Jaccard similarity > jaccard_threshold.
        - Expand with high-prevalence mutations (overlap_threshold).
        - If resulting clade has >= min_clade_size, return it.
        - Otherwise, return None.

        Returns
        -------
        (mutations_list, child_row_indices) or None
        """
        # Enumerate mutations and their sets
        all_muts = _enumerate_mutations(parent_idx, truncal_muts)
        if not all_muts:
            return None, None

        mut_sets: dict[Mutation, set[int]] = {}
        mut_sets = {m: _mutation_set(parent_idx, m) for m in all_muts}
        sorted_muts = sorted(
            (m for m in mut_sets if len(mut_sets[m]) > min_clade_size),
            key=lambda m: len(mut_sets[m]),
            reverse=True,
        )

        for seed in sorted_muts:
            base_set = mut_sets[seed]

            # Build initial split mutation set
            split_muts = {seed}
            for m in sorted_muts:
                if _jaccard(base_set, mut_sets[m]) > jaccard_threshold:
                    split_muts.add(m)
                # if len(split_muts) >= min_mutations:
                #    break

            if len(split_muts) < min_mutations:
                continue

            # Expand split mutation set with high-prevalence mutations
            child_idx = _match_set_with_missing(parent_idx, split_muts, max_mismatch_frac=0)
            split_muts = _high_prevalence_mutations(child_idx, threshold=overlap_threshold)
            # split_muts.update(additional_muts)
            child_idx = _match_set_with_missing(parent_idx, split_muts, max_mismatch_frac=0.02)

            if len(child_idx) >= min_clade_size:
                return split_muts, child_idx

        return None, None

    def _split_clade(
        node_id: str,
        parent_idx: np.ndarray,
        truncal_muts: set[Mutation],
        depth: int = 0,
    ) -> None:
        """
        Recursively split this clade:

        - While we can find valid splits (>= min_clade_size), peel off children.
        - Each split produces a child that is further split (with truncal muts excluded).
        - After no more splits mark remaining as doublet if < doublet_threshold of total,
        """
        remaining = set(parent_idx)
        n_splits = 0

        while len(remaining) >= min_clade_size:
            current_idx = np.fromiter(remaining, dtype=int)
            clade_min_mutations = min_mutations if (n_splits == 0 or depth == 0) else min_mutations - 1
            split_muts, child_idx = _find_best_split(current_idx, truncal_muts, clade_min_mutations)
            if (split_muts is None) or ((len(child_idx) == len(parent_idx)) & depth > 0):
                break

            # Peel off child clade
            child_cell_ids = cell_ids[child_idx]
            child_id = next(node_name_gen)
            tree.add_node(child_id, descendants=list(child_cell_ids), depth=depth + 1, stump=True)
            tree.add_edge(node_id, child_id)

            # Recurse on child with updated used set
            if depth < max_depth - 1:
                _split_clade(child_id, child_idx, split_muts, depth + 1)
            n_splits += 1

            # Remove those cells from remaining pool
            remaining -= set(child_idx)

        # Mark remaining as doublet if small enough
        if (len(remaining) < len(parent_idx) * doublet_threshold) and remaining:
            residual_idx = np.fromiter(remaining, dtype=int)
            child_cell_ids = cell_ids[residual_idx]
            child_id = next(node_name_gen)
            tree.add_node(child_id, descendants=list(child_cell_ids), doublet=True, depth=depth + 1, stump=True)
            tree.add_edge(node_id, child_id)
        # Otherwise, stop splitting
        elif len(remaining) > 0:
            tree.remove_nodes_from(nx.descendants(tree, node_id))

        # Stop splitting if >3 splits detected
        if depth > 0 and n_splits >= 3:
            tree.remove_nodes_from(nx.descendants(tree, node_id))

    # Root clade contains all cells
    root_id = next(node_name_gen)
    tree.add_node(root_id, descendants=list(cell_ids), depth=0, stump=True)
    _split_clade(root_id, np.arange(n_cell, dtype=int), truncal_muts=set(), depth=0)
    doublets = np.array(
        [
            cell_id
            for node in tree.nodes
            if tree.nodes[node].get("doublet", False)
            for cell_id in tree.nodes[node]["descendants"]
        ]
    )
    if connect_leaves:
        tree = _connect_leaves(tree)
        for node in tree.nodes:
            if "descendants" in tree.nodes[node]:
                del tree.nodes[node]["descendants"]

    return tree, doublets


def fasttree(
    characters: pd.DataFrame,
    mutation_rate: float = 0.1,
    root_name: str = "root",
    root_method: str | None = "character_outgroup",
    missing_state: str = "-",
    node_name_gen: callable = node_name_generator(),
    fasttree_cmd: str = "FastTree",
    logfile: str | None = None,
) -> nx.DiGraph:
    """Reconstruct a lineage tree with FastTree.

    Parameters
    ----------
    characters : pd.DataFrame
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
    number_of_states = len(set(characters.values.ravel().tolist()) - {missing_state})
    if number_of_states > 20:
        raise ValueError("number of unique states must be <= 20.")

    with tempfile.TemporaryDirectory() as td:
        fasta_path = os.path.join(td, "character_matrix.fasta")
        edit_base = os.path.join(td, "edit_distance")
        tree_path = os.path.join(td, "tree.nwck")

        save_edit_distance(number_of_states, edit_base)
        save_characters_fasta(characters, fasta_path, missing_state=missing_state)

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
            tree = newick_to_tree(newick_str, node_name_gen=node_name_gen)
        # Root tree
        if root_method is None:
            pass
        elif root_method == "centroid":
            tree = centroid(tree, new_root=root_name)
        elif root_method == "character_outgroup":
            tree = character_outgroup(tree, characters, new_root=root_name)
        else:
            raise ValueError(f"Unknown root_method: {root_method}")

    return tree


def _resolve_polytomies(
    tree: nx.DiGraph,
    node_name_gen: callable = node_name_generator(),
) -> nx.DiGraph:
    """Arbitrarily resolve polytomies so every internal node has <= 2 children.

    LAML-Pro requires a rooted *binary* tree as input. This peels off children
    of any node with > 2 children into a caterpillar of new internal nodes.
    """
    tree = tree.copy()
    for node in list(tree.nodes):
        children = list(tree.successors(node))
        while len(children) > 2:
            c1 = children.pop()
            c2 = children.pop()
            new_node = next(node_name_gen)
            tree.add_node(new_node)
            for c in (c1, c2):
                attrs = tree.get_edge_data(node, c) or {}
                tree.remove_edge(node, c)
                tree.add_edge(new_node, c, **attrs)
            tree.add_edge(node, new_node)
            children.append(new_node)
    return tree


def laml(
    characters: pd.DataFrame,
    initial_tree: nx.DiGraph,
    mode: str = "search",
    mutation_priors: dict[str, float] | dict[str, dict[str, float]] | None = None,
    ultrametric: bool = True,
    max_iterations: int = 2500,
    temp: float = 0.1,
    min_branch_length: float = 0.01,
    seed: int = 73,
    threads: int = 10,
    missing_state: str = "-",
    unedited_state: str = "*",
    node_name_gen: callable = node_name_generator(),
    lamlpro_cmd: str = "lamlpro",
    logfile: str | None = None,
) -> nx.DiGraph:
    """Reconstruct a lineage tree with LAML-Pro (maximum likelihood PMMO model).

    LAML-Pro requires an initial rooted tree (e.g. from :func:`fasttree`). Its
    topology and branch lengths are then optimized (``mode="search"``) or only the
    branch lengths are fit (``mode="optimize"``) under the PMMO model. The
    LAML-Pro Python API is not used (it is buggy); the ``lamlpro`` command-line
    tool is invoked instead.

    Parameters
    ----------
    characters : pd.DataFrame
        String-encoded character matrix. Index = cell IDs, columns = sites.
        ``unedited_state`` = no mutation, ``missing_state`` = missing, anything
        else = a specific mutation state.
    initial_tree : nx.DiGraph
        Initial rooted tree. Its leaves must match ``characters.index`` (order does
        not matter). Polytomies are resolved arbitrarily, since LAML-Pro requires a
        binary tree.
    mode : str
        ``"search"`` to optimize topology + branch lengths, ``"optimize"`` to fit
        only branch lengths for the initial topology.
    mutation_priors : dict, optional
        Prior mutation probabilities, either as ``{state: probability}`` (applied to
        every character) or ``{column: {state: probability}}`` (per character).
        States must use the same string encoding as ``characters``.
    ultrametric : bool
        Enforce equal root-to-leaf path lengths.
    max_iterations : int
        Maximum hill-climbing iterations (search mode).
    temp : float
        Starting temperature for the simulated-annealing topology search.
    min_branch_length : float
        Minimum branch length relative to a scaled tree of unit height.
    seed : int
        Random seed for reproducibility.
    threads : int
        Number of threads.
    missing_state, unedited_state : str
        Encoding of missing and unedited states.
    node_name_gen : Callable
        Generator yielding unique node IDs.
    lamlpro_cmd : str
        Path to / name of the ``lamlpro`` binary.
    logfile : Optional[str]
        If set, console output from LAML-Pro is appended to this file.

    Returns
    -------
    nx.DiGraph
        Rooted tree with branch lengths stored on the ``"length"`` edge attribute
        and cumulative node times on the ``"time"`` node attribute. Each LAML-Pro
        branch length is the time before that node, so the root's time is its own
        stem branch length and each child's time adds the branch leading into it.
    """
    if mode not in {"search", "optimize"}:
        raise ValueError("mode must be 'search' or 'optimize'.")
    number_of_states = len(set(characters.values.ravel().tolist()) - {missing_state, unedited_state})
    if number_of_states == 0:
        raise ValueError("characters contains no mutations.")

    # Validate that the initial tree's leaves match the character matrix.
    leaves = set(get_leaves(initial_tree))
    cells = set(characters.index)
    if leaves != cells:
        raise ValueError(
            "initial_tree leaves must match characters.index "
            f"({len(leaves - cells)} extra leaves, {len(cells - leaves)} missing)."
        )

    # LAML-Pro requires a rooted binary tree: collapse unifurcations (out_degree 1)
    # and resolve polytomies (out_degree > 2) so every internal node bifurcates.
    initial_tree = collapse_unifurcations(initial_tree.copy(), allow_root=True)
    initial_tree = _resolve_polytomies(initial_tree, node_name_gen=node_name_gen)
    input_root = get_root(initial_tree)

    with tempfile.TemporaryDirectory() as td:
        matrix_path = os.path.join(td, "character_matrix.csv")
        tree_path = os.path.join(td, "initial_tree.nwk")
        priors_path = os.path.join(td, "mutation_priors.csv")
        out_prefix = os.path.join(td, "lamlpro")
        out_tree_path = f"{out_prefix}_tree.newick"

        mapping = save_characters_csv(
            characters, matrix_path, missing_state=missing_state, unedited_state=unedited_state
        )
        # Write internal node names so that, in optimize mode, LAML-Pro preserves
        # them and the output node names match the input tree.
        with open(tree_path, "w") as f:
            f.write(tree_to_newick(initial_tree, record_node_names=True) + "\n")

        priors_flag = ""
        if mutation_priors is not None:
            # A flat {state: prob} dict is applied to every character; a nested
            # {column: {state: prob}} dict specifies priors per character.
            if all(not isinstance(v, dict) for v in mutation_priors.values()):
                mutation_priors = dict.fromkeys(characters.columns, mutation_priors)
            col_to_idx = {col: i for i, col in enumerate(characters.columns)}
            with open(priors_path, "w") as f:
                for col, states in mutation_priors.items():
                    for state, prob in states.items():
                        f.write(f"{col_to_idx[col]},{mapping[state]},{prob}\n")
            priors_flag = f"--mutation-priors {priors_path} "

        cmd = (
            f"bash -c '. ~/.bashrc && "
            f"{lamlpro_cmd} "
            f"--matrix {matrix_path} "
            f"--tree {tree_path} "
            f"--output {out_prefix} "
            f"--data-type character-matrix "
            f"--mode {mode} "
            f"--max-iterations {max_iterations} "
            f"--temp {temp} "
            f"--min-branch-length {min_branch_length} "
            f"--seed {seed} "
            f"--threads {threads} "
            f"{priors_flag}"
            f"{'--ultrametric ' if ultrametric else ''}'"
        )

        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if logfile:
            with open(logfile, "ab") as lf:
                lf.write(stdout or b"")
                lf.write(stderr or b"")
        if process.returncode != 0:
            raise RuntimeError(f"LAML-Pro failed (code {process.returncode}). Stderr:\n{stderr.decode('utf-8')}")

        with open(out_tree_path) as f:
            newick_str = f.read().strip()
        tree = newick_to_tree(newick_str, node_name_gen=node_name_gen)
        # LAML-Pro assigns each node a branch length equal to the time before that
        # node. The root's branch length (the stem before the first division) is not
        # stored as an edge, so read it directly to seed the root's time.
        root_length = float(ete3.Tree(newick_str, format=1).dist or 0.0)

    # LAML-Pro preserves the input rooting; tidy it up and restore the root name.
    tree = collapse_unifurcations(tree, allow_root=True)
    root = get_root(tree)
    if root != input_root:
        tree = nx.relabel_nodes(tree, {root: input_root})

    # Record branch lengths as cumulative node "time" attributes. The root's time is
    # its own stem branch length; each child adds the branch length leading into it.
    times = {input_root: root_length}
    for parent, child in nx.bfs_edges(tree, input_root):
        times[child] = times[parent] + tree[parent][child].get("length", 0.0)
    nx.set_node_attributes(tree, times, "time")

    return tree
