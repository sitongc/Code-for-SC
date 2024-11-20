#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statistics as st
import warnings
import sys
import os
import random
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -i [input] -o [output file] -s [sample id]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-i", "--input", help="path of h5ad input")
    optparser.add_option("-o", "--out", help="output file")
    optparser.add_option("-s", "--sample", help="sample id")
    optparser.add_option("-m", "--meta", help="meta info")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.meta or not options.out:
        optparser.print_help()
        sys.exit(-1)
        
    print("loading the sample...")
    #adata = sc.read_visium(options.input)
    adata = sc.read_h5ad(options.input)
    adata.var_names_make_unique()
    
    meta = pd.read_csv(options.meta)
    
    print("filtering the cell with extremely low expression...")
    filter_barcode = meta["barcode"]
    adata = adata[filter_barcode].copy()

    print("clustering...")

    random.seed(100)

    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="clusters", resolution = 1)
    print("Done")
    
    adata.obs["SingleR_orig"] = meta["singleR_orig"].values
    #adata.obs["majority_vote"] = meta["majority_vote"].values

    print("Generating a UMAP plot...")
    output = [options.out, options.sample, "_UMAP.pdf"]
    output_path = "/".join(str(i) for i in output)
    plt.rcParams["figure.figsize"] = (4, 4)
    with plt.rc_context():  # Use this to set figure params like size and dpi
        sc.pl.umap(adata, color=["SingleR_orig","clusters"], wspace=0.4, show = False)
        plt.savefig(output_path, bbox_inches="tight")
    
    print("Generating a spatial plot...")
    output = [options.out, options.sample, "_spatial.pdf"]
    output_path = "/".join(str(i) for i in output)
    plt.rcParams["figure.figsize"] = (6, 4)
    with plt.rc_context():  # Use this to set figure params like size and dpi
        sc.pl.spatial(adata, img_key="hires", color=["SingleR_orig","clusters"], size=1.5, show = False)
        plt.savefig(output_path, bbox_inches="tight")
        
    output = [options.out, options.sample, "data.h5ad"]
    output_path = "/".join(str(i) for i in output)
    
    adata.write(output_path)       

    #Quality control
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    #pre Quality control
    output = [options.out, options.sample, "pre_QC.pdf"]
    output_path = "/".join(str(i) for i in output)
    plt.rcParams["figure.figsize"] = (6, 4)
    with plt.rc_context():  # Use this to set figure params like size and dpi
         fig, axs = plt.subplots(1, 3, figsize=(15, 4))
         sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
         sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1])
         sns.histplot(adata.obs["pct_counts_mt"], kde=False, bins=60, ax=axs[2])
         plt.savefig(output_path, bbox_inches="tight")


    #after Quality control 
    data_sd = st.stdev(adata.obs["total_counts"])
    data_mean = st.mean(adata.obs["total_counts"])

    lower_cutoff = data_mean - 2*data_sd
    upper_cutoff = data_mean + 2*data_sd

    print(f"#Spot before filtering: {adata.n_obs}")
    sc.pp.filter_cells(adata, min_counts=lower_cutoff)
    sc.pp.filter_cells(adata, max_counts=upper_cutoff)
    adata = adata[adata.obs["pct_counts_mt"] < 10].copy()

    print(f"#Spot after filtering: {adata.n_obs}")
    sc.pp.filter_genes(adata, min_cells=10)
    
    output = [options.out, options.sample, "after_QC.pdf"]
    output_path = "/".join(str(i) for i in output)
    plt.rcParams["figure.figsize"] = (6, 4)
    with plt.rc_context():  # Use this to set figure params like size and dpi
         fig, axs = plt.subplots(1, 3, figsize=(15, 4))
         sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
         sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1])
         sns.histplot(adata.obs["pct_counts_mt"], kde=False, bins=60, ax=axs[2])
         plt.savefig(output_path, bbox_inches="tight")



    #save the cluster info and UMAP info
    output = [options.out, options.sample, "cluster_leiden.csv"]
    output_path = "/".join(str(i) for i in output)
    cluster_info = pd.DataFrame(adata.obs["clusters"], index=adata.obs_names)
    cluster_info.to_csv(output_path)

    umap_coord = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
    output = [options.out, options.sample, "umap_coord.csv"]
    output_path = "/".join(str(i) for i in output)
    umap_coord.to_csv(output_path)

    output = [options.out, options.sample, "data_processed.h5ad"]
    output_path = "/".join(str(i) for i in output)

    adata.write(output_path)
    print("Finish")


if __name__=='__main__':
    main()
