import os
import subprocess
import tempfile

import networkx as nx
import numpy as np
import pandas as pd

from .config import node_name_generator
from .root import centroid, character_outgroup
from .utils import newick_to_tree, save_characters_fasta, save_edit_distance

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
            if (split_muts is None) or (len(child_idx) == len(parent_idx)):
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
        else:
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
