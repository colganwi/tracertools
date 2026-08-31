import multiprocessing as mp
import re
from collections.abc import Sequence

import Levenshtein
import numpy as np
import pandas as pd
import pysam
import tqdm as tqdm
from sklearn.mixture import GaussianMixture

from .config import edit_ids


def insertion_from_alignment(sequence, cigar, pos, ref_begin=0, window=2):
    """Extract insertion sequence from a read alignment."""
    if pos < ref_begin:
        return None
    pos = pos - ref_begin
    ref_pos = 0
    seq_pos = 0
    insertions = ""
    for op in re.findall(r"(\d+)([MIDNSHP])", cigar):
        length, operation = int(op[0]), op[1]
        if operation == "M":
            seq_pos += length
            ref_pos += length
        elif operation == "S":
            seq_pos += length
        elif operation == "D":
            ref_pos += length
        elif operation == "I":
            if (pos - window <= ref_pos) and (ref_pos <= pos + window):
                indel_pos = seq_pos - ref_pos + pos
                insertions += sequence[indel_pos : indel_pos + length]
            seq_pos += length
        if ref_pos > pos + window:
            break
    if ref_pos < pos:
        return None
    elif insertions == "":
        return "*"
    return insertions


def barcode_from_alignment(sequence, cigar, start, stop, ref_begin=0):
    """Extract barcode sequence from a read alignment."""
    aligned_sequence = ""
    start = start - ref_begin
    stop = stop - ref_begin
    seq_pos = 0
    ref_pos = 0
    for op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length, operation = int(op[0]), op[1]
        if operation in ["M", "S"]:
            if ref_pos + length > start:
                sub_start = max(0, start - ref_pos)
                sub_end = min(length, stop - ref_pos)
                aligned_sequence += sequence[seq_pos + sub_start : seq_pos + sub_end]
            seq_pos += length
            if operation == "M":
                ref_pos += length
        elif operation == "D":
            ref_pos += length
        elif operation == "I":
            seq_pos += length
    return aligned_sequence


def correct_to_whitelist(seqs, whitelist, max_dist=2):
    """Correct sequences to a whitelist using Levenshtein distance"""
    unique_seqs = list(dict.fromkeys(seqs))
    corrections = {}
    for seq in unique_seqs:
        closest_match = None
        min_distance = float("inf")
        for white_seq in whitelist:
            distance = Levenshtein.distance(seq, white_seq)
            if distance <= max_dist and distance < min_distance:
                closest_match = white_seq
                min_distance = distance
        corrections[seq] = closest_match if closest_match is not None else seq
    corrected_seqs = [corrections[str] for str in seqs]
    return corrected_seqs


def sigma_threshold(x, max_sigma=3, n_components=2, min_threshold=1, log=False):
    """Estimate threshold based on Gaussian Mixture Model."""
    # If more than 5000 values sample
    if len(x) > 5000:
        x = x.sample(5000)
    if log:
        x = np.log1p(x)
    gmm = GaussianMixture(n_components=n_components, random_state=0).fit(x.values.reshape(-1, 1))
    high_component = np.argmax(gmm.means_)
    threshold = gmm.means_[high_component, 0] - max_sigma * np.sqrt(gmm.covariances_[high_component, 0, 0])
    if log:
        threshold = np.expm1(threshold)
    return max(min_threshold, threshold)


def select_allele(allele, sites=["RNF2", "HEK3", "EMX1"]):
    """Select allele given conflicting sequencing reads."""
    agg_funcs = {col: "sum" if col in ["UMI", "readCount", "frac"] else "first" for col in allele.columns}
    aggregated = allele.groupby("n_alleles").agg(agg_funcs)
    n_edits = 0
    for site in sites:
        values = list(allele[site].dropna().unique())
        if len(values) == 1:
            continue
        elif len(values) == 2 and "None" in values:
            aggregated[site] = values[0] if values[1] == "None" else values[1]
            n_edits += 1
        elif len(values) > 1:
            return allele
    if n_edits == 1:
        aggregated["n_alleles"] = 1
        return aggregated
    else:
        return allele


def resolve_alleles(df: pd.DataFrame, *, sites: Sequence[str] | None = None) -> pd.DataFrame:
    """Resolve alleles with conflicting sequencing reads (memory optimized)."""
    if df.empty:
        return df

    if sites is None:
        default_sites = ["RNF2", "HEK3", "EMX1"]
        sites = [c for c in default_sites if c in df.columns]
    else:
        sites = [c for c in sites if c in df.columns]

    keys = ["intID", "cellBC"]
    if not sites:
        # still need to add n_alleles like original
        n_alleles = df.groupby(keys, sort=False)["intID"].transform("size")
        out = df.copy(deep=False)
        out["n_alleles"] = n_alleles
        return out

    # Compute n_alleles WITHOUT copying whole df first
    n_alleles = df.groupby(keys, sort=False)["intID"].transform("size")
    out = df.copy(deep=False)  # shallow view; we only add one column
    out["n_alleles"] = n_alleles

    # Work only on the n_alleles==2 subset
    g2 = out[out["n_alleles"] == 2]
    if g2.empty:
        return out

    # Grouped nunique on just the sites columns (smaller than whole frame)
    g2g = g2.groupby(keys, sort=False)

    # total uniques (including '*')
    u_total = g2g[sites].nunique(dropna=True)

    # uniques excluding '*' (treat '*' as NA) computed without copying all of g2
    # This creates an intermediate sites-only block, not a full df copy.
    sites_block = g2[sites].mask(g2[sites].eq("*"))
    u_non_none = sites_block.groupby([g2[k] for k in keys], sort=False).nunique(dropna=True)
    u_non_none.index.names = keys

    # whether the pair contains literal '*' for each site (vectorized)
    has_none = g2g[sites].agg(lambda s: (s == "*").any())

    # resolvable if pair is {x, '*'} at that site
    resolvable = (u_total.eq(2)) & (u_non_none.eq(1)) & has_none

    # any conflict anywhere (>1 distinct non-'*' values)
    conflict_any = (u_non_none > 1).any(axis=1)

    # good pairs: exactly one resolvable site, no conflicts elsewhere
    edit_counts = resolvable.sum(axis=1)
    good_pairs = (~conflict_any) & (edit_counts == 1)
    if not bool(good_pairs.any()):
        return out

    good_idx = resolvable.index[good_pairs]  # MultiIndex of (intID, cellBC)

    # Identify rows in g2 that belong to the "good" groups, without building a big temp df
    g2_idx = pd.MultiIndex.from_frame(g2[keys], names=keys)
    idx_mask = g2_idx.isin(good_idx)

    to_collapse = g2.loc[idx_mask]
    unresolved_pairs = g2.loc[~idx_mask]
    not_two = out[out["n_alleles"] != 2]

    # Aggregate other columns
    num_cols = [c for c in ("UMI", "readCount", "frac") if c in to_collapse.columns]
    excluded = set(keys) | set(sites) | {"n_alleles"} | set(num_cols)
    keep_first_cols = [c for c in to_collapse.columns if c not in excluded]

    agg_dict = dict.fromkeys(num_cols, "sum")
    for c in keep_first_cols:
        agg_dict[c] = "first"

    collapsed = to_collapse.groupby(keys, sort=False).agg(agg_dict)

    # Pick site values with minimal extra memory:
    # mask '*' -> NA, take first non-NA per group, fill remaining with '*'
    for s in sites:
        picked = to_collapse[s].mask(to_collapse[s] == "*").groupby([to_collapse[k] for k in keys], sort=False).first()
        picked.index.names = keys
        collapsed[s] = picked.reindex(collapsed.index).fillna("*").to_numpy()

    collapsed["n_alleles"] = 1
    collapsed = collapsed.reset_index()

    # Final concat
    return pd.concat([not_two, unresolved_pairs, collapsed], ignore_index=True, copy=False)


def _init_worker(df, sites):
    global _DF, _SITES
    _DF = df
    _SITES = sites


def _resolve_alleles_worker(intID):
    # No copy here; resolve_alleles only shallow-copies when adding a column.
    sub = _DF[_DF["intID"] == intID]
    return resolve_alleles(sub, sites=_SITES)


def resolve_alleles_parallel(
    df: pd.DataFrame,
    sites=None,
    threads: int = 10,
    chunksize: int = 20,
    concat_batch: int = 50,
) -> pd.DataFrame:
    """Resolve alleles in parallel with minimized peak memory."""
    unique_ids = df["intID"].unique()

    # Prefer fork when available (Linux) to leverage copy-on-write and avoid pickling df.
    # On macOS/Windows this will fall back to spawn semantics.
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        ctx = mp.get_context()

    with ctx.Pool(
        processes=threads,
        initializer=_init_worker,
        initargs=(df, sites),
        maxtasksperchild=200,  # helps long runs avoid fragmentation
    ) as pool:
        it = pool.imap(_resolve_alleles_worker, unique_ids, chunksize=chunksize)

        # Chunked concat: reduces peak RAM vs list(results) then single concat
        out_chunks = []
        buf = []
        for res in tqdm.tqdm(it, total=len(unique_ids), desc="Resolving alleles"):
            buf.append(res)
            if len(buf) >= concat_batch:
                out_chunks.append(pd.concat(buf, ignore_index=True, copy=False))
                buf.clear()

        if buf:
            out_chunks.append(pd.concat(buf, ignore_index=True, copy=False))

    return pd.concat(out_chunks, ignore_index=True, copy=False)


def read_sam(file, verbose=False):
    """Read SAM file into a DataFrame."""
    samfile = pysam.AlignmentFile(file, "rb")
    records = []
    for alignment in samfile.fetch(until_eof=True):
        record_data = {
            "query_name": alignment.query_name,
            "is_read1": alignment.is_read1,
            "is_read2": alignment.is_read2,
            "ref": alignment.reference_name if alignment.reference_name else "unmapped",
            "query_begin": alignment.query_alignment_start,
            "ref_begin": alignment.reference_start,
            "map_qual": alignment.mapping_quality,
            "CIGAR": alignment.cigarstring,
            "seq": alignment.query_sequence,
        }
        records.append(record_data)
    samfile.close()
    df = pd.DataFrame(records)
    n = len(df)
    n_unmapped = sum(df["ref"] == "unmapped")
    if verbose:
        print(f"{n} reads with {n_unmapped} unmapped ({n_unmapped / n * 100:.2f}%)")
    return df


def alleles_to_characters(alleles, edit_ids=edit_ids, min_prob=None, other_id="!", order=None, index="cellBC"):
    """Convert allele table to character matrix"""
    characters = alleles.copy()
    if isinstance(index, str):
        index = [index]
    # Map alleles to characters
    for site, mapping in edit_ids.items():
        characters[site] = characters[site].map(mapping).fillna(other_id)
        if min_prob is not None and f"{site}_prob" in characters.columns:
            characters.loc[characters[f"{site}_prob"] < min_prob, site] = "-"
    characters = pd.melt(
        characters[index + ["intID"] + list(edit_ids.keys())],
        id_vars=index + ["intID"],
        var_name="site",
        value_name="allele",
    )
    characters = characters.pivot_table(
        index=index, columns=["intID", "site"], values="allele", aggfunc="first"
    ).fillna("-")

    # sort by max allele fraction
    def max_fraction(int_id):
        int_data = characters.xs(int_id, level=0, axis=1)
        counts = int_data.apply(pd.Series.value_counts, axis=0).fillna(0)
        total_counts = counts.sum(axis=0)
        valid_counts = counts.loc[lambda x: x.index > "-"]  # Exclude - and *
        max_fraction_value = (valid_counts / total_counts).max().max()
        return max_fraction_value

    if order is None:
        order = sorted(characters.columns.levels[0], key=max_fraction, reverse=True)
    characters = characters.reindex(order, level=0, axis=1)

    # Reindex
    characters.columns = [f"{intID}-{site}" for intID, site in characters.columns]
    return characters
