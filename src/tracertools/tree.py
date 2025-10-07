import os
import subprocess
import tempfile
from typing import Optional

import cassiopeia as cas
import networkx as nx
import numpy as np
import pandas as pd
import pycea
import treedata as td

from .utils import bfs_names, mask_truncal_edits, newick_to_tree, save_characters_fasta, save_edit_distance


def get_root(tree):
    """Finds the root of a tree"""
    return [node for node in tree.nodes if tree.in_degree(node) == 0][0]


def get_leaves(tree: nx.DiGraph):
    """Finds the leaves of a tree"""
    return [node for node in nx.dfs_postorder_nodes(tree, get_root(tree)) if tree.out_degree(node) == 0]


def infer_ancestral_characters(tdata,tree = "tree",key = "characters",edit_cost = .6,copy = False):
    """Reconstruct ancestral characters using Sankoff algorithm"""
    tree_key = tree
    if copy:
        tdata = tdata.copy()
    n_characters = max(tdata.obsm[key].max().max() + 1,10)
    costs = np.ones(shape = (n_characters,n_characters),dtype=float)
    costs[0,:] = edit_cost
    np.fill_diagonal(costs,0)
    costs = pd.DataFrame(costs,index = range(0,n_characters),columns = range(0,n_characters))
    pycea.tl.ancestral_states(tdata,keys = key,method = "sankoff",costs = costs,missing_state=-1,tree = tree_key)
    if copy:
        return tdata


def collapse_mutationless_edges(tdata,tree = "tree",key = "characters",tree_added = "tree",mutation_key = None,copy = False):
    """Collapse mutationless edges"""
    tree_key = tree
    if copy:
        tdata = tdata.copy()
    tree = tdata.obst[tree_key].copy()
    root = [node for node in tree.nodes if tree.in_degree(node) == 0][0]
    for edge in reversed(list(nx.dfs_edges(tree,root))):
        if mutation_key is not None:
            has_mutation = tree.edges[edge][mutation_key]
        else:
            has_mutation = np.any(tree.nodes[edge[1]][key] != tree.nodes[edge[0]][key])
        if not has_mutation:
            children = list(tree.successors(edge[1]))
            if len(children) > 0:
                for child in children:
                    tree.add_edge(edge[0],child)
                tree.remove_edge(*edge)
                tree.remove_node(edge[1])
    tdata.obst[tree_added] = tree
    if copy:
        return tdata


def majority_character(characters,min_size = 20,min_frac = .8):
    """Find the majority character in a list of characters"""
    characters = np.array(characters)
    if len(characters) < min_size:
        return -1
    characters = characters[characters != -1]
    if len(characters) == 0:
        return -1
    unique_values, counts = np.unique(characters, return_counts=True)
    for value, count in zip(unique_values, counts):
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


def identify_mutations(tdata,tree = "tree",key = "characters",key_added = "has_mutation",min_frac = .75,copy = False):
    """Mark edges with a mutation"""
    if copy:
        tdata = tdata.copy()
    method = lambda x: majority_character(x,min_frac = min_frac)
    pycea.tl.ancestral_states(tdata, keys = key, method = method,keys_added="majority_characters",tree = tree)
    for edge in tdata.obst[tree].edges:
        has_mutation = not (same_characters(tdata.obst[tree].nodes[edge[0]]["majority_characters"],
                                          tdata.obst[tree].nodes[edge[1]]["majority_characters"]) and
                            same_characters(tdata.obst[tree].nodes[edge[0]][key],
                                            tdata.obst[tree].nodes[edge[1]][key]))
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
    logfile: Optional[str] = None,
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
        raise ValueError(
            f"Detected {len(detected)} states, but number_of_states={number_of_states}."
        )

    with tempfile.TemporaryDirectory() as td:
        fasta_path = os.path.join(td, "character_matrix.fasta")
        edit_base = os.path.join(td, "edit_distance")
        tree_path = os.path.join(td, "tree.nwck")

        save_edit_distance(number_of_states, edit_base)
        save_characters_fasta(character_matrix, fasta_path)

        cmd = (
            f". ~/.bashrc && "
            f"{fasttree_cmd} -nosupport -noml -nome "
            f"-matrix {edit_base} "
            f"{fasta_path} > {tree_path}"
        )

        process = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if logfile:
            with open(logfile, "ab") as lf:
                lf.write(stdout or b"")
                lf.write(stderr or b"")
        if process.returncode != 0:
            raise RuntimeError(
                f"FastTree failed (code {process.returncode}). Stderr:\n{stderr.decode('utf-8')}"
            )
        # Read and (optionally) root with ete3
        with open(tree_path, "r") as f:
            newick_str = f.read().strip()
            tree = newick_to_tree(newick_str, midpoint_root=True if add_root else False)

    return tree


def reconstruct_cassiopeia(characters,solver = "upgma"):
    """Reconstruct a tree using cassiopeia"""
    solvers = {"nj":cas.solver.NeighborJoiningSolver(cas.solver.dissimilarity.weighted_hamming_distance,
                                                 add_root=True,fast = False),
           "upgma":cas.solver.UPGMASolver(cas.solver.dissimilarity.weighted_hamming_distance,fast = False),
           "greedy":cas.solver.VanillaGreedySolver()}
    cas_tree = cas.data.CassiopeiaTree(character_matrix=characters)
    solvers[solver].solve(cas_tree)
    tree = cas_tree.get_tree_topology()
    tree = nx.relabel_nodes(tree,bfs_names(tree))
    return tree


def reconstruct_tree(tdata,solver = "fasttree",key = "characters",tree_added = "tree",
                     mask_truncal = True, logfile = None, copy = False):
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


def estimate_leaf_fitness(tdata,tree = "tree",depth_key = "depth",key_added = "fitness",copy = False):
    """Estimate leaf fitness based on tree and branch lengths"""
    tree_key = tree
    if copy:
        tdata = tdata.copy()
    nx_tree = tdata.obst[tree_key].copy()
    for node in nx_tree:
        nx_tree.nodes[node]["_depth"] = nx_tree.nodes[node][depth_key]
    cas_tree = cas.data.CassiopeiaTree(tree = nx_tree)
    for edge in cas_tree.depth_first_traverse_edges():
        t1 = cas_tree.get_attribute(edge[0],"_depth")
        t2 = cas_tree.get_attribute(edge[1],"_depth")
        cas_tree.set_branch_length(edge[0], edge[1], abs(t1-t2))
    fitness_estimator = cas.tools.fitness_estimator.LBIJungle()
    fitness_estimator.estimate_fitness(cas_tree)
    fitnesses = np.array([cas_tree.get_attribute(cell, 'fitness') for cell in cas_tree.leaves])
    fitnesses = pd.Series(fitnesses, index=cas_tree.leaves)
    tdata.obs[key_added] = tdata.obs_names.map(fitnesses)
    if copy:
        return tdata



def estimate_branch_lengths(tdata: td.TreeData,tree: str = "tree",key: str = "characters",minimum_branch_length = .001,
                              key_added: str = "time",solver: str = "CLARABEL",verbose: bool = False,copy: bool = False):
    """Estimate branch lengths using IIDExponentialMLE

    Parameters
    ----------
    tdata : the TreeData object
    tree : the key for the tree in tdata.obst
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
    tree_key = tree
    if copy:
        tdata = tdata.copy()
    t = tdata.obst[tree_key]
    cas_tree = cas.data.CassiopeiaTree(character_matrix = tdata.obsm[key],tree = t)
    cas_tree.set_all_character_states(nx.get_node_attributes(cas_tree.get_tree_topology(),key))
    IIDExponentialMLE(verbose=verbose, solver = solver, minimum_branch_length=minimum_branch_length).estimate_branch_lengths(cas_tree)
    node_times = nx.get_node_attributes(cas_tree.get_tree_topology(),"time")
    nx.set_node_attributes(t,node_times,key_added)
    del tdata.obst[tree]
    tdata.obst[tree] = t
    if copy:
        return tdata
