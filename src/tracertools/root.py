
import networkx as nx
import pandas as pd

from .tree import collapse_unifurcations
from .utils import get_leaves, get_root


def reroot_on_edge(G: nx.DiGraph, new_root, best_edge):
    """Insert a new root on best_edge and orient edges away from it."""
    U = G.to_undirected()
    u, v = best_edge

    # Work on a copy of the undirected tree
    U2 = U.copy()
    if U2.has_edge(u, v):
        U2.remove_edge(u, v)
    U2.add_node(new_root)
    U2.add_edge(new_root, u)
    U2.add_edge(new_root, v)

    # Build a directed tree oriented away from new_root
    H = nx.DiGraph()
    H.add_nodes_from(G.nodes(data=True))
    H.add_node(new_root)

    for p, c in nx.bfs_edges(U2, new_root):
        # copy attrs if present in either direction in original G
        attrs = (G.get_edge_data(p, c) or G.get_edge_data(c, p) or {})
        H.add_edge(p, c, **attrs)

    # Remove old root & any other unifurcations
    H = collapse_unifurcations(H, allow_root=True)
    return H


def centroid(G: nx.DiGraph, new_root="root"):
    """Reroot a tree to minimize leaf imbalance.

    1) Find the edge whose two sides have the closest leaf counts.
    2) Insert a new root on that edge and orient edges away from it.
    3) Collapse unifurcations (out_degree == 1), including the old root.
    Returns: (H, chosen_edge, (L1, L2), new_root)
    """
    U = G.to_undirected()
    if not nx.is_tree(U):
        raise ValueError("Input must be a tree.")

    # --- Compute leaf counts via a single DFS ---
    r = next(iter(U.nodes))
    parent = {r: None}
    order = []
    stack = [r]
    while stack:
        u = stack.pop()
        order.append(u)
        for v in U[u]:
            if v == parent.get(u):
                continue
            parent[v] = u
            stack.append(v)

    is_leaf = {u: (U.degree[u] == 1) for u in U}
    leaves_down = {u: 1 if is_leaf[u] else 0 for u in U}
    for u in reversed(order):
        for v in U[u]:
            if v == parent.get(u):
                continue
            leaves_down[u] += leaves_down[v]
    total_leaves = sum(is_leaf.values())

    def edge_leaf_sides(u, v):
        # Return (Lu, Lv) = leaf counts on the two sides of edge {u,v}
        if parent.get(v) == u:
            Lv = leaves_down[v]
            Lu = total_leaves - Lv
        elif parent.get(u) == v:
            Lu = leaves_down[u]
            Lv = total_leaves - Lu
        else:  # shouldn’t happen in a tree, but keep symmetric fallback
            Lu, Lv = leaves_down[u], leaves_down[v]
        return Lu, Lv

    # --- Pick best edge: minimize |Lu-Lv|; tie-break by maximizing min(Lu,Lv), deterministic id ---
    best = None
    best_edge = None
    for u, v in U.edges():
        Lu, Lv = edge_leaf_sides(u, v)
        diff = abs(Lu - Lv)
        key = (diff, -min(Lu, Lv), tuple(sorted((u, v), key=str)))
        if best is None or key < best:
            best, best_edge = key, (u, v)

    # --- Reroot using helper ---
    return reroot_on_edge(G, new_root, best_edge)


def _get_character_outgroup(characters: pd.DataFrame, missing_state = "-", unedited_state = "*"):
    """Identify outgroup based on characters."""
    best_col = None
    best_value = None
    best_prop = -1
    for col in characters.columns:
        col_values = characters[col][~characters[col].isin([missing_state, unedited_state])]
        if col_values.empty:
            continue
        mode_value = col_values.mode().iloc[0]
        prop = (characters[col] == mode_value).mean()
        if prop > best_prop:
            best_prop = prop
            best_col = col
            best_value = mode_value
    if best_col is None:
        return []
    return characters.index[characters[best_col] == best_value].tolist()


def character_outgroup(tree: nx.DiGraph, characters, new_root="root"):
    """Reroot tree based on character outgroup."""
    t = tree.copy()
    outgroup = _get_character_outgroup(characters)
    leaves = get_leaves(t)
    root = get_root(t)
    total_leaves = len(leaves)
    total_outgroup = len(outgroup)
    nx.set_node_attributes(t, {leaf: {"total": 1,"outgroup": 0} for leaf in leaves})
    nx.set_node_attributes(t, {leaf: {"outgroup": 1} for leaf in outgroup})
    for node in nx.dfs_postorder_nodes(t, source=root):
        if node in leaves:
            continue
        children = list(t.successors(node))
        total = sum(t.nodes[child]["total"] for child in children)
        outgroup = sum(t.nodes[child]["outgroup"] for child in children)
        t.nodes[node]["total"] = total
        t.nodes[node]["outgroup"] = outgroup
    def _node_jaccard(node):
        return max(t.nodes[node]["outgroup"] / (t.nodes[node]["total"] + total_outgroup - t.nodes[node]["outgroup"]),
                     (total_outgroup - t.nodes[node]["outgroup"]) /
                     (total_leaves - t.nodes[node]["total"] + total_outgroup - (total_outgroup - t.nodes[node]["outgroup"])))
    best_jaccard = -1
    best_edge = None
    for u, v in t.edges():
        child_jaccard = _node_jaccard(v)
        parent_jaccard = _node_jaccard(u)
        jaccard = max(child_jaccard, parent_jaccard)
        if jaccard > best_jaccard or (jaccard == best_jaccard and parent_jaccard < child_jaccard):
            best_jaccard = jaccard
            best_edge = (u, v)
    return reroot_on_edge(tree, new_root, best_edge)
