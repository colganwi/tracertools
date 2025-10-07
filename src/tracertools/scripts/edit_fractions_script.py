import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def get_parser():
    """Get parser for edit_fractions script"""
    parser = argparse.ArgumentParser(description="Calculate edit fractions from allele counts files", add_help=False)
    # Add arguments
    parser.add_argument("--path", type=str, default="./data/", help="Path to directory allele counts")
    parser.add_argument("--sites", type=str, default="EMX1,HEK3,RNF2", help="Site to calculate edit fractions for")
    parser.set_defaults(func=edit_fractions)
    return parser


def edit_fractions(path, sites, **kwargs):
    """Calculate edit fractions from allele counts files"""
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
    # Calculate edit fractions
    edit_counts = allele_counts[sites + ["readCount", "sample", "intID"]].melt(
        id_vars=["readCount", "sample", "intID"], value_name="edit", var_name="site"
    )
    edit_frac = (
        edit_counts.query("edit != '-'").groupby(["sample", "site"])["readCount"].sum()
        / edit_counts.groupby(["sample", "site"])["readCount"].sum()
    ).reset_index(name="edit_frac")
    edit_frac.to_csv(base_path / "edit_fractions.csv")
    # Plot edit fractions
    site_palette = {"EMX1": "#1874CD", "HEK3": "#CD2626", "RNF2": "#FFE600"}
    sns.barplot(data=edit_frac, x="sample", y="edit_frac", hue="site", palette=site_palette, saturation=1)
    plt.xticks(rotation=90)
    plt.savefig(base_path / "edit_fractions.png", bbox_inches="tight")
    plt.close()
