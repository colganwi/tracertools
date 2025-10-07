import re

import Levenshtein
import numpy as np
import pandas as pd
import pysam
from sklearn.mixture import GaussianMixture
from typing import Sequence

from .config import edit_ids


def insertion_from_alignment(sequence, cigar, pos, ref_begin = 0, window=2):
    """Extract insertion sequence from a read alignment."""
    if pos < ref_begin:
        return None
    pos = pos - ref_begin
    ref_pos = 0
    seq_pos = 0
    insertions = ""
    for op in re.findall(r'(\d+)([MIDNSHP])', cigar):
        length, operation = int(op[0]), op[1]
        if operation == 'M':
            seq_pos += length
            ref_pos += length
        elif operation == "S":
            seq_pos += length
        elif operation == 'D':
            ref_pos += length
        elif operation == 'I':
            if (pos - window <= ref_pos) and (ref_pos <= pos + window):
                indel_pos = seq_pos - ref_pos + pos
                insertions += sequence[indel_pos:indel_pos + length]
            seq_pos += length
        if ref_pos > pos + window:
            break
    if ref_pos < pos:
        return None
    elif insertions == "":
        return "-"
    return insertions

def barcode_from_alignment(sequence, cigar, start, stop, ref_begin = 0):
    """Extract barcode sequence from a read alignment."""
    aligned_sequence = ""
    start = start - ref_begin
    stop = stop - ref_begin
    seq_pos = 0
    ref_pos = 0
    for op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length, operation = int(op[0]), op[1]
        if operation in ["M","S"]:
            if ref_pos + length > start:
                sub_start = max(0, start - ref_pos)
                sub_end = min(length, stop - ref_pos)
                aligned_sequence += sequence[seq_pos + sub_start:seq_pos + sub_end]
            seq_pos += length
            if operation == "M":
                ref_pos += length
        elif operation == "D":
            ref_pos += length
        elif operation == "I":
            seq_pos += length
    return aligned_sequence


def correct_to_whitelist(seqs, whitelist, max_dist = 2):
    """Correct sequences to a whitelist using Levenshtein distance"""
    unique_seqs = list(dict.fromkeys(seqs))
    corrections = {}
    for seq in unique_seqs:
        closest_match = None
        min_distance = float('inf')
        for white_seq in whitelist:
            distance = Levenshtein.distance(seq, white_seq)
            if distance <= max_dist and distance < min_distance:
                closest_match = white_seq
                min_distance = distance
        corrections[seq] = closest_match if closest_match is not None else seq
    corrected_seqs = [corrections[str] for str in seqs]
    return corrected_seqs


def sigma_threshold(x, max_sigma = 3, n_components = 2, min_threshold = 1,log = False):
    """Estimate threshold based on Gaussian Mixture Model."""
    # If more than 5000 values sample
    if len(x) > 5000:
        x = x.sample(5000)
    if log:
        x = np.log1p(x)
    gmm = GaussianMixture(n_components=n_components, random_state=0).fit(x.values.reshape(-1, 1)) 
    high_component = np.argmax(gmm.means_)
    threshold = gmm.means_[high_component,0] - max_sigma * np.sqrt(gmm.covariances_[high_component,0,0])
    if log:
        threshold = np.expm1(threshold)
    return max(min_threshold, threshold)


def select_allele(allele, sites=["RNF2", "HEK3", "EMX1"]):
    """Select allele given conflicting sequencing reads."""
    agg_funcs = {col: 'sum' if col in ["UMI", "readCount", "frac"] else 'first' for col in allele.columns}
    aggregated = allele.groupby("n_alleles").agg(agg_funcs)
    n_edits = 0
    for site in sites:
        values = list(allele[site].dropna().unique())
        if len(values) == 1:
            continue
        elif len(values) == 2 and 'None' in values:
            aggregated[site] = values[0] if values[1] == 'None' else values[1]
            n_edits += 1
        elif len(values) > 1:
            return allele
    if n_edits == 1:
        aggregated["n_alleles"] = 1
        return aggregated
    else:
        return allele


def resolve_alleles(df: pd.DataFrame, *, sites: Sequence[str] | None = None) -> pd.DataFrame:
    """Resolve alleles with conflicting sequencing reads."""
    if sites is None:
        default_sites = ["RNF2", "HEK3", "EMX1"]
        sites = [c for c in default_sites if c in df.columns]
    keys = ["intID", "cellBC"]
    out = df.copy()
    gc_keys = out.groupby(keys, sort=False)
    out["n_alleles"] = gc_keys["intID"].transform("size")
    g2 = out[out["n_alleles"] == 2].copy()
    if g2.empty or not sites:
        return out
    g2_gc = g2.groupby(keys, sort=False)
    # total uniques
    u_total = g2_gc[sites].nunique(dropna=True)
    g2_masked = g2.copy()
    for s in sites:
        g2_masked[s] = g2_masked[s].where(g2_masked[s] != "-")
    u_non_none = g2_masked.groupby(keys, sort=False)[sites].nunique(dropna=True)
    # does the pair contain a literal '-' for each site?
    has_none = g2.groupby(keys, sort=False)[sites].apply(lambda df_: df_.eq("-").any(axis=0))
    # resolvable if pair is {x, '-'}
    resolvable = (u_total.eq(2)) & (u_non_none.eq(1)) & has_none
    # any conflict (>1 distinct non-'None' values at any site)?
    conflict_any = (u_non_none > 1).any(axis=1)
    # good pairs: exactly one resolvable site, no conflicts elsewhere
    edit_counts = resolvable.sum(axis=1)
    good_pairs = (~conflict_any) & (edit_counts == 1)
    if not good_pairs.any():
        return out
    good_idx = resolvable.index[good_pairs]
    idx_mask = g2.set_index(keys).index.isin(good_idx)
    to_collapse = g2.loc[idx_mask]
    num_cols = [c for c in ["UMI", "readCount", "frac"] if c in g2.columns]
    excluded = set(keys + list(sites) + ["n_alleles"] + num_cols)
    keep_first_cols = [c for c in g2.columns if c not in excluded]
    agg_dict = {c: "sum" for c in num_cols}
    agg_dict.update({c: "first" for c in keep_first_cols})
    collapsed = (
        to_collapse.groupby(keys, sort=False)
        .agg(agg_dict)
        .reset_index()
    )

    def pick_site(col: str) -> pd.Series:
        subset = to_collapse[[*keys, col]]
        pref = (
            subset[subset[col] != "-"]
            .dropna()
            .drop_duplicates(subset=keys)
            .set_index(keys)[col]
        )
        uniq = (
            subset.dropna()
            .drop_duplicates(subset=keys + [col])
            .groupby(keys, sort=False)[col].first()
        )
        out_series = uniq.copy()
        out_series.loc[pref.index] = pref
        return out_series

    for s in sites:
        collapsed[s] = pick_site(s).reindex(collapsed.set_index(keys).index).values
    collapsed["n_alleles"] = 1
    unresolved_pairs = g2.loc[~idx_mask]
    not_two = out[out["n_alleles"] != 2]
    result = pd.concat([not_two, unresolved_pairs, collapsed], ignore_index=True)
    return result

def read_sam(file,verbose=False):
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
        print(f"{n} reads with {n_unmapped} unmapped ({n_unmapped/n*100:.2f}%)")
    return df


def alleles_to_characters(alleles,edit_ids = edit_ids,min_prob = None,other_id = 9,order = None,index = "cellBC"):
    """Convert allele table to character matrix"""
    characters = alleles.copy()
    if isinstance(index,str):
        index = [index]
    # Map alleles to characters
    for site, mapping in edit_ids.items():
        characters[site] = characters[site].map(mapping).fillna(other_id).astype(int)
        if min_prob is not None and f"{site}_prob" in characters.columns:
            characters.loc[characters[f"{site}_prob"] < min_prob,site] = -1
    characters = pd.melt(characters[index + ["intID"] + list(edit_ids.keys())],
                               id_vars = index + ["intID"],var_name = "site",value_name = "allele")
    characters = characters.pivot_table(index = index,columns = ["intID","site"],values = "allele").fillna(-1).astype(int)
    # sort by max allele fraction
    def max_fraction(int_id):
        int_data = characters.xs(int_id, level=0, axis=1)
        counts = int_data.apply(pd.Series.value_counts, axis=0).fillna(0)
        total_counts = counts.sum(axis=0)
        valid_counts = counts.loc[lambda x: x.index > 0]  # Exclude -1 and 0
        max_fraction_value = (valid_counts / total_counts).max().max()
        return max_fraction_value
    if order is None:
        order = sorted(characters.columns.levels[0], key=max_fraction, reverse=True)
    characters = characters.reindex(order, level=0, axis=1)
    # Reindex
    characters.columns = ['{}-{}'.format(intID, site) for intID, site in characters.columns]
    return characters
