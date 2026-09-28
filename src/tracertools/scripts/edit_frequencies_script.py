import argparse
from pathlib import Path

import pandas as pd


def get_parser():
    """Get parser for edit_frequencies script"""
    parser = argparse.ArgumentParser(description="Calculate edit frequencies from allele counts files", add_help=False)
    # Add arguments
    parser.add_argument("--path", type=str, default="./data/", help="Path to directory allele counts")
    parser.add_argument("--sites", type=str, default="EMX1,HEK3,RNF2", help="Site to calculate edit frequencies for")
    parser.add_argument("--top_n", type=int, default=8, help="Number of top edits to keep")
    parser.set_defaults(func=edit_frequencies)
    return parser


def edit_frequencies(path, sites, top_n, **kwargs):
    """Calculate edit frequencies from allele counts files"""
    path = Path(path)
    base_path = path.parent
    sites = sites.split(",")
    # Load data
    allele_counts = []
    for file in list(path.glob("*.csv")):
        if file.name.startswith("."):
            continue
        sample = Path(file).stem.replace("_allele_counts", "")
        allele_counts.append(pd.read_csv(file).assign(sample=sample))
    allele_counts = pd.concat(allele_counts, ignore_index=True)
    # Calculate edit frequencies
    long_counts = allele_counts.melt(
        id_vars=[c for c in allele_counts.columns if c not in ["RNF2", "HEK3", "EMX1"]],
        value_vars=["RNF2", "HEK3", "EMX1"],
        var_name="site",
        value_name="edit",
    )
    edit_counts = long_counts.groupby(["site", "edit", "sample"])["readCount"].sum().reset_index()
    edit_counts["rank"] = edit_counts.groupby(["site", "sample"])["readCount"].rank(method="first", ascending=False)
    edit_counts["edit"] = edit_counts.apply(lambda row: row["edit"] if row["rank"] <= top_n else "other", axis=1)
    edit_counts = edit_counts.groupby(["site", "edit", "sample"])["readCount"].sum().reset_index()
    edit_counts.to_csv(base_path / "edit_frequencies.csv")
