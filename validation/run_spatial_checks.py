"""Synthetic integration checks; run from github-preparation/."""
import ast
import importlib.metadata
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

root = Path.cwd()
work = Path(tempfile.mkdtemp(prefix="spatial_code_checks_", dir="/private/tmp"))
package = work / "carotid-cta-remodelling"
shutil.copytree(root / "carotid-cta-remodelling", package)
(package / "data").mkdir(exist_ok=True)
gene_dir = package / "results/02_concordant_gene_sets"
gene_dir.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(2026)
x = np.log1p(rng.poisson(3, (120, 1200)).astype(float))
adata = ad.AnnData(x, obs=pd.DataFrame({"region": ["A", "B"] * 60,
                                      "cell_type": ["synthetic"] * 120},
                                     index=[f"cell{i}" for i in range(120)]),
                   var=pd.DataFrame(index=[f"gene{i}" for i in range(1200)]))
# A different raw matrix catches accidental default scoring from adata.raw.
adata.raw = ad.AnnData(rng.random(x.shape), obs=adata.obs.copy(), var=adata.var.copy())
adata.write_h5ad(package / "data/public_xenium_panel2.h5ad")
(gene_dir / "concordant_positive_RI_symbols.txt").write_text("gene1\ngene2\ngene1\n")
(gene_dir / "concordant_negative_RI_symbols.txt").write_text("missing_gene\n")
script = package / "scripts/06_public_xenium_spatial_score.py"
ast.parse(script.read_text())
result = subprocess.run([os.sys.executable, str(script)], cwd=package, text=True, capture_output=True)
if result.returncode:
    raise RuntimeError(result.stderr)
output = package / "results/06_public_xenium_spatial_score"
actual = pd.read_csv(output / "public_xenium_RI_module_scores_per_cell.csv")
expected = adata.copy()
sc.tl.score_genes(expected, ["gene1", "gene2"], score_name="expected", use_raw=False, random_state=2026)
assert np.allclose(actual.positive_RI_spatial_score, expected.obs.expected)
assert actual.negative_RI_spatial_score.isna().all()
assert (output / "public_xenium_RI_module_scores_by_region.csv").exists()
assert (output / "public_xenium_RI_module_scores_by_celltype.csv").exists()
adata.X = rng.poisson(3, x.shape).astype(float)
adata.write_h5ad(package / "data/public_xenium_panel2.h5ad")
bad = subprocess.run([os.sys.executable, str(script)], cwd=package, text=True, capture_output=True)
assert bad.returncode != 0 and "raw integer counts" in bad.stderr
print("PASS spatial scoring against direct Scanpy reference, explicit X, duplicate genes, empty gene set")
print("PASS raw-count rejection")
print("Scanpy", importlib.metadata.version("scanpy"))
print("Synthetic fixtures:", work)
