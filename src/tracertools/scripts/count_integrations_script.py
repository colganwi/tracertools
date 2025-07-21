import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from tracertools.seq import sigma_threshold


def get_parser():
    """Get parser for count_integrations script"""
    parser = argparse.ArgumentParser(description="Count integrations from allele counts files", add_help=False)
    # Add arguments
    parser.add_argument("--path", type=str, default="./data/", help="Path to directory allele counts")
    parser.add_argument("--threshold", type=int, default=None, help="Read count threshold for filtering integrations")
    parser.set_defaults(func=count_integrations)
    return parser


def count_integrations(path, threshold=None, **kwargs):
    """Count integrations from allele counts files"""
    # Set paths
    path = Path(path)
    base_path = path.parent
    # Load data
    allele_counts = []
    for file in list(path.glob("*.csv")):
        if file.name.startswith("."):
            continue
        sample = Path(file).stem.replace("_allele_counts", "")
        allele_counts.append(pd.read_csv(file).assign(sample=sample))
    allele_counts = pd.concat(allele_counts, ignore_index=True)
    # Set threshold
    if threshold is None:
        threshold = sigma_threshold(allele_counts["readCount"][allele_counts["readCount"] > 10], log=True)
    fig, ax = plt.subplots(figsize=(3, 3))
    sns.histplot(allele_counts, x="readCount", log_scale=True, bins=30, ax=ax)
    ax.set_yscale("log")
    ax.axvline(threshold, color="black", linestyle="--")
    fig.savefig(base_path / "read_count_threshold.png", bbox_inches="tight")
    # Filter counts
    filtered_counts = allele_counts[allele_counts["readCount"] > threshold]
    # Calculate number of integrations per sample
    n_int = filtered_counts.groupby("sample").agg(n=("intID", "nunique"), intIDs=("intID", lambda x: ",".join(x)))
    n_int.to_csv(base_path / "integrations_per_sample.csv")
    # Plot integrations per sample
    counts_wide = filtered_counts.pivot(index="intID", columns="sample", values="readCount").fillna(0)
    counts_wide = np.log10(counts_wide + 1)
    sns.clustermap(counts_wide)
    plt.savefig(base_path / "integration_clustermap.png", bbox_inches="tight")
