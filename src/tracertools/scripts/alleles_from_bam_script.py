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
    parser = argparse.ArgumentParser(description="Call alleles from bam file",add_help=False)
    # Add arguments
    parser.add_argument("--bam", type=str,help="Bam file")
    parser.add_argument("--out", type=str, help="Output csv file")
    parser.add_argument("--barcode_start", type=int, default=270, help="Start of integration barcode in reference")
    parser.add_argument("--barcode_end", type=int, default=300, help="End of integration barcode in reference")
    parser.add_argument("--site_positions", type=str, default="{'RNF2':332,'HEK3':380,'EMX1':448}", help="Dictionary mapping site names to positions")
    parser.add_argument("--min_reads", type=int, default=2, help="Minimum number of reads to call allele")
    parser.add_argument("--extract_barcode", type = bool, default=False, help="Extract barcode sequence from reads instead of using alignment")
    parser.set_defaults(func=alleles_from_bam)
    return parser

def has_cb_tag(bamfile):
    i = 0
    for read in bamfile.fetch(until_eof=True):
        i += 1
        if i > 1000:
            break
        if read.has_tag("CB"):
            return True
    return False


def call_alleles_10x(param):
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


def call_alleles_bulk(bam, barcode_start, barcode_end, sites, min_reads=10, extract_barcode=False):
    # Process reads
    bamfile = pysam.AlignmentFile(bam, "rb")
    molecules = defaultdict(lambda: { **{ name: None for name in sites.keys() }, "intID": None})
    total_reads = bamfile.count()
    bamfile.reset()
    for read in tqdm(bamfile.fetch(until_eof=True), total=total_reads, mininterval=20, desc="TS"):
        if (read.is_secondary or read.mapping_quality < 30):
            continue
        if read.reference_start < barcode_start and read.reference_end > barcode_start:
            molecules[read.query_name]["intID"] = read.reference_name
        for name, pos in sites.items():
            allele = insertion_from_alignment(
                    read.query_sequence,
                    read.cigarstring,
                    pos,
                    read.reference_start
                )
            if allele is not None:
                molecules[read.query_name][name] = allele
    # Aggregate alleles
    alleles = pd.DataFrame.from_dict(molecules, orient='index')
    alleles = alleles.dropna(axis=0, how='any')
    for col in sites.keys():
        alleles = alleles[~alleles[col].str.contains("N")]
    if "EMX1" in alleles.columns:
        alleles["EMX1"] = alleles["EMX1"].replace("CTTGGG", "-")
    alleles["readCount"] = 1
    alleles = alleles.groupby(["intID"] + list(sites.keys())).agg({"readCount":"sum"}).reset_index()
    alleles = alleles.query(f"readCount >= {min_reads}")
    return alleles


def alleles_from_bam(bam,
                    out,
                    barcode_start,
                    barcode_end,
                    site_positions,
                    min_reads=10,
                    extract_barcode=False,
                    **kwargs):
     # Parse the arguments
    sites = ast.literal_eval(site_positions)
    bamfile = pysam.AlignmentFile(bam, "rb")
    intIDs = [ref for ref in bamfile.references if "intID" in ref]
    has_cb = has_cb_tag(bamfile)
    bamfile.close()
    # All alleles 10x
    if has_cb:
        lock = mp.Manager().Lock()
        pd.DataFrame(columns=["intID", "cellBC"] + list(sites.keys()) + ["UMI", "readCount"]).to_csv(out, index=False)
        if extract_barcode:
            call_alleles_10x((intIDs[0],bamfile,barcode_start,barcode_end,sites,min_reads,extract_barcode,out,lock))
        else:
            with mp.Pool(processes=8) as pool:
                _ = list(tqdm(pool.imap_unordered(call_alleles_10x,
                [(intID,bamfile,barcode_start,barcode_end,sites,min_reads,extract_barcode,out,lock) for intID in intIDs]),
                total=len(intIDs),mininterval=60, desc="TS"))
    # Call alleles bulk
    else:
        alleles = call_alleles_bulk(bam, barcode_start, barcode_end, sites, min_reads, extract_barcode)
        alleles.to_csv(out, index=False)
