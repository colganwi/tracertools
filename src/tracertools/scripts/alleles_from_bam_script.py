import argparse
import ast
import multiprocessing as mp
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam
from tqdm.auto import tqdm
from tracertools.seq import barcode_from_alignment, insertion_from_alignment


def get_parser():
    """Get parser for align_experiments script"""
    parser = argparse.ArgumentParser(description="Call alleles from bam file")
    # Add arguments
    parser.add_argument("--bam", type=str,help="Bam file")
    parser.add_argument("--out", type=str, help="Output csv file")
    parser.add_argument("--barcode_start", type=int, help="Start of integration barcode in reference")
    parser.add_argument("--barcode_end", type=int, help="End of integration barcode in reference")
    parser.add_argument("--site_positions", type=str, help="Dictionary mapping site names to positions")
    parser.add_argument("--min_reads", type=int, help="Minimum number of reads to call allele")
    parser.add_argument("--extract_barcode", type = bool, help="Extract barcode sequence from reads instead of using alignment")
    parser.set_defaults(func=call_alleles)
    return parser

def call_alleles(param):
    # Setup
    intID, bam, barcode_start, barcode_end, sites, min_reads, extract_barcode, out, lock = param
    if len(sites) > 0:
        end = max(sites.values())
    else:
        end = barcode_end
    # Get iterator
    bamfile = pysam.AlignmentFile(bam, "rb")
    if extract_barcode:
        total_reads = bamfile.count(contig=intID)
        read_iter = tqdm(bamfile.fetch(intID), total=total_reads, mininterval=60, desc="TS")
    else:
        read_iter = bamfile.fetch(intID)
    # Process reads
    umi_counts = defaultdict(int)
    for read in read_iter:
        if (read.mapping_quality < 30 or
            read.reference_start > barcode_start or
            read.reference_end < end or
            'N' in read.cigarstring or
            not read.has_tag('CB') or
            not read.has_tag('UB')):
            continue
        # Get integration
        if extract_barcode:
            intID = barcode_from_alignment(read.query_sequence, read.cigarstring, barcode_start, barcode_end, read.reference_start)
        else:
            intID = read.reference_name
        # Get allele
        alleles = []
        for name, pos in sites.items():
            allele = insertion_from_alignment(read.query_sequence, read.cigarstring, pos, read.reference_start)
            if name == "EMX1" and allele == "CTTGGG":
                allele = "None"
            alleles.append(allele)
        key = (read.get_tag('UB'),read.get_tag('CB'), intID, *alleles)
        umi_counts[key] += 1
    bamfile.close()
    # Correct and aggregate UMIs
    if len(umi_counts) > 0:
        site_names = list(sites.keys())
        allele_counts = pd.DataFrame(umi_counts.keys(), columns=["UMI","cellBC","intID"] + site_names)
        allele_counts["readCount"] = umi_counts.values()
        del umi_counts
        # correct UMIs
        allele_counts = allele_counts.groupby(["intID","cellBC","UMI"] + site_names).agg(
            {"readCount":"sum"}).sort_values("readCount", ascending=False).reset_index()
        agg_dict = {site: 'first' for site in site_names} 
        agg_dict["readCount"] = "sum"
        allele_counts = allele_counts.groupby(["intID","cellBC","UMI"]).agg(agg_dict).reset_index()
        # collapse UMIs
        allele_counts = allele_counts.groupby(["intID","cellBC"] + site_names).agg(
        {"UMI":"size","readCount":"sum"}).reset_index()
        # filter alleles
        allele_counts = allele_counts.query(f"readCount >= {min_reads}")
        with lock:
            allele_counts.to_csv(out, mode='a', header=False, index=False)

def alleles_from_bam(bam,
                    out,
                    barcode_start,
                    barcode_end,
                    site_positions,
                    min_reads=10,
                    extract_barcode=False):
     # Parse the arguments
    sites = ast.literal_eval(site_positions)
    bamfile = pysam.AlignmentFile(bam, "rb")
    intIDs = [ref for ref in bamfile.references if "intID" in ref]
    bamfile.close()
    # Make output file
    lock = mp.Manager().Lock()
    pd.DataFrame(columns=["intID", "cellBC"] + list(sites.keys()) + ["UMI", "readCount"]).to_csv(out, index=False)
    # Process
    if extract_barcode:
        call_alleles((intIDs[0],bamfile,barcode_start,barcode_end,sites,min_reads,extract_barcode,out,lock))
    # Process in parallel
    else:
        with mp.Pool(processes=8) as pool:
            _ = list(tqdm(pool.imap_unordered(call_alleles,
            [(intID,bamfile,barcode_start,barcode_end,sites,min_reads,extract_barcode,out,lock) for intID in intIDs]),
            total=len(intIDs),mininterval=60, desc="TS"))
