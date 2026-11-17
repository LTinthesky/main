#!/usr/bin/env python3

"""
Public Xenium spatial module score analysis.
Requires an AnnData .h5ad file with region annotations.
"""

from pathlib import Path

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

# -------------------------
# User parameters
# -------------------------
h5ad_file = Path("data/public_xenium_panel2.h5ad")
positive_gene_file = Path("results/02_concordant_gene_sets/concordant_positive_RI_symbols.txt")
negative_gene_file = Path("results/02_concordant_gene_sets/concordant_negative_RI_symbols.txt")
output_dir = Path("results/06_public_xenium_spatial_score")

region_col = "region"
celltype_col = "cell_type"
positive_score_col = "positive_RI_spatial_score"
negative_score_col = "negative_RI_spatial_score"

output_dir.mkdir(parents=True, exist_ok=True)

# -------------------------
# Load data and gene sets
# -------------------------
adata = sc.read_h5ad(h5ad_file)
if not adata.var_names.is_unique or not adata.obs_names.is_unique:
    raise ValueError("Gene and cell IDs must be unique")
# Explicitly use X, not the optional adata.raw matrix. X must be log-normalized
# across a suitable background gene pool before this script is run.
values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
if not np.isfinite(values).all():
    raise ValueError("Expression X must contain finite values")
if values.size and np.allclose(values, np.round(values)):
    raise ValueError("X appears to contain raw integer counts; provide log-normalized expression")

positive_genes = [x.strip() for x in positive_gene_file.read_text().splitlines() if x.strip()]
negative_genes = [x.strip() for x in negative_gene_file.read_text().splitlines() if x.strip()]

positive_genes = sorted(set(positive_genes).intersection(adata.var_names))
negative_genes = sorted(set(negative_genes).intersection(adata.var_names))

(output_dir / "mapped_positive_RI_spatial_genes.txt").write_text("\n".join(positive_genes) + "\n")
(output_dir / "mapped_negative_RI_spatial_genes.txt").write_text("\n".join(negative_genes) + "\n")

# -------------------------
# Scanpy score_genes
# -------------------------
if len(positive_genes) > 0:
    sc.tl.score_genes(adata, gene_list=positive_genes, score_name=positive_score_col, use_raw=False, random_state=2026)
else:
    adata.obs[positive_score_col] = np.nan

if len(negative_genes) > 0:
    sc.tl.score_genes(adata, gene_list=negative_genes, score_name=negative_score_col, use_raw=False, random_state=2026)
else:
    adata.obs[negative_score_col] = np.nan

# -------------------------
# Export summaries
# -------------------------
obs_cols = [positive_score_col, negative_score_col]
if region_col in adata.obs.columns:
    obs_cols.append(region_col)
if celltype_col in adata.obs.columns:
    obs_cols.append(celltype_col)

score_df = adata.obs[obs_cols].copy()
score_df.insert(0, "cell_id", adata.obs_names)
score_df.to_csv(output_dir / "public_xenium_RI_module_scores_per_cell.csv", index=False)

if region_col in score_df.columns:
    region_summary = (
        score_df
        .groupby(region_col, dropna=False)
        .agg(
            n_cells=(positive_score_col, "size"),
            positive_median=(positive_score_col, "median"),
            positive_mean=(positive_score_col, "mean"),
            negative_median=(negative_score_col, "median"),
            negative_mean=(negative_score_col, "mean"),
        )
        .reset_index()
        .sort_values("positive_median", ascending=False)
    )
    region_summary.to_csv(output_dir / "public_xenium_RI_module_scores_by_region.csv", index=False)

if celltype_col in score_df.columns:
    celltype_summary = (
        score_df
        .groupby(celltype_col, dropna=False)
        .agg(
            n_cells=(positive_score_col, "size"),
            positive_median=(positive_score_col, "median"),
            positive_mean=(positive_score_col, "mean"),
            negative_median=(negative_score_col, "median"),
            negative_mean=(negative_score_col, "mean"),
        )
        .reset_index()
        .sort_values("positive_median", ascending=False)
    )
    celltype_summary.to_csv(output_dir / "public_xenium_RI_module_scores_by_celltype.csv", index=False)

print(f"Saved public Xenium spatial score results to: {output_dir}")
