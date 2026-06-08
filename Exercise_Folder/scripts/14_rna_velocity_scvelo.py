#!/usr/bin/env python3
"""
Module 14 (Python) — RNA velocity with scVelo.

Companion to the R script 14_trajectory_cellcomm.R (Slingshot pseudotime +
CellChat communication) and the laptop notebook
Exercise_Folder/Tutorial_14_Trajectory_CellCommunication.qmd.

RNA velocity is a Python analysis (scVelo), so it lives in this .py file rather
than the .R script. It needs a SPLICED/UNSPLICED count matrix — a .loom from
`velocyto run` or `kb count --workflow lamanno`. The standard Cell Ranger
*filtered* matrix does NOT contain the unspliced layer, so velocity cannot be
run from it (the most common reason a velocity run "fails" or returns flat
arrows).

Run (e.g. on Talapas, inside a Python env that has scvelo installed):
    OUT_DIR=../output/Mod14 LOOM=../data/ifnb_velocyto.loom python 14_rna_velocity_scvelo.py

Install once (choose one):
    pip install scvelo scanpy
    mamba install -c conda-forge -c bioconda scvelo scanpy
"""

import os

import scanpy as sc
import scvelo as scv

OUT_DIR = os.environ.get("OUT_DIR", "../output/Mod14")
LOOM = os.environ.get("LOOM", "../data/ifnb_velocyto.loom")
os.makedirs(OUT_DIR, exist_ok=True)

scv.settings.verbosity = 1
scv.settings.figdir = OUT_DIR  # scv.pl.*(save=...) writes into this directory

# 1. Load spliced/unspliced counts (from velocyto / kb — NOT cellranger filtered).
adata = scv.read(LOOM)

# 2. Standard scVelo preprocessing.
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# 3. Make sure a UMAP embedding exists to project velocities onto. If your loom
#    already carries the Seurat UMAP (transferred via obsm["X_umap"]), this is
#    skipped; otherwise compute a quick one here.
if "X_umap" not in adata.obsm:
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    sc.tl.umap(adata)

# 4. Estimate velocities with the dynamical model (Bergen et al. 2020).
scv.tl.recover_dynamics(adata, n_jobs=4)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# 5. Figure out: velocity stream on the UMAP embedding -> ../output/Mod14/.
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    title="RNA velocity (scVelo, dynamical) — Module 14",
    save="Mod14_velocity_stream.png",
    show=False,
)

# 6. Save the AnnData (with velocity layers) for any downstream use.
adata.write(os.path.join(OUT_DIR, "Mod14_scvelo.h5ad"))

print(f"Wrote Mod14 RNA-velocity figure + h5ad to {OUT_DIR}")
