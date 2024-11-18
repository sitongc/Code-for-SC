#!/usr/bin/env python

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statistics as st
import celltypist

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

#input the data 
adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
#adata = sc.read_h5ad('V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()

# QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# filter the cells based on the cutoff
data_sd = st.stdev(adata.obs["total_counts"])
data_mean = st.mean(adata.obs["total_counts"])

lower_cutoff = data_mean - 2*data_sd
upper_cutoff = data_mean + 2*data_sd

sc.pp.filter_cells(adata, min_counts=lower_cutoff)
sc.pp.filter_cells(adata, max_counts=upper_cutoff)
adata = adata[adata.obs["pct_counts_mt"] < 10].copy()

# clustering
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters", resolution = 1)

# annotation by celltypist
model = celltypist.models.Model.load("Immune_All_Low.pkl")
predictions = celltypist.annotate(adata, model=model, majority_voting=True)
adata.obs["cell_type"] = predictions.predicted_labels["predicted_labels"]

# visualization 

#UMAP
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
sc.pl.umap(adata, color="cell_type", title="Human Lymph", size=80, legend_fontsize=10)

# spatial
plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
sc.pl.spatial(adata, img_key="hires",color="cell_type", legend_fontsize=10)



