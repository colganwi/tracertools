import re

import Levenshtein
import numpy as np
import pandas as pd
import pysam
from sklearn.mixture import GaussianMixture


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
