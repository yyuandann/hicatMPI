import pickle
import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import re
import ast
import time

from transcriptomic_clustering.final_merging import final_merge, FinalMergeKwargs

def main(adata_path, latent_path, out_dir, final_merge_kwargs): # data is a tuple of anndata and key argumnents
    start = time.perf_counter()

    def preprocess_dict(dict_str):
    # Remove comments using regex
        dict_str_cleaned = re.sub(r"#.*", "", dict_str)
        return dict_str_cleaned
    
    final_merge_kwargs = preprocess_dict(final_merge_kwargs)
    final_merge_kwargs = ast.literal_eval(final_merge_kwargs)

    pca_kwargs = final_merge_kwargs['pca_kwargs']
    filter_pcs_kwargs = final_merge_kwargs['filter_pcs_kwargs']
    filter_known_modes_kwargs = final_merge_kwargs['filter_known_modes_kwargs']
    project_kwargs = final_merge_kwargs['project_kwargs']
    merge_clusters_kwargs = final_merge_kwargs['merge_clusters_kwargs']
    latent_kwargs = final_merge_kwargs['latent_kwargs']
    
    merge_kwargs = FinalMergeKwargs(
        pca_kwargs = pca_kwargs,
        filter_pcs_kwargs = filter_pcs_kwargs,
        filter_known_modes_kwargs = filter_known_modes_kwargs,
        project_kwargs = project_kwargs,
        merge_clusters_kwargs = merge_clusters_kwargs,
        latent_kwargs = latent_kwargs
    )

    if '.zarr' in adata_path:
        print(f"Reading in anndata from zarr format: {adata_path}", flush=True)
        adata = ad.read_zarr(adata_path)
    elif '.h5ad' in adata_path:
        adata = sc.read(adata_path)
    else:
        raise ValueError(f"Unsupported file format for adata_path: {adata_path}. Please provide a .h5ad or .zarr file.")
    
    print(f"Finished reading in anndata: {adata}", flush=True)

    if np.max(adata.X) > 100:
        print(f"Raw count data provided")
        print(f"Normlazing total counts to 1e6...")
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        print(f"Finished normalization. max:{np.max(adata.X)}")
    else:
        print(f"Normalized data provided")

    # determine if latent_path is valid
    if os.path.exists(latent_path):
        latent = pd.read_csv(latent_path, index_col=0)
        latent.index = latent.index.astype(str)
        latent = latent.loc[adata.obs_names]
        adata.obsm['latent'] = np.asarray(latent)
        merge_kwargs.latent_kwargs['latent_component'] = 'latent'
        print(f"Finished reading in latent space: {latent.shape}")
    else:
        print(f"scvi_path is invalid or not provided. Using latent_kwargs['latent_component'] for clustering (None for PCA and a str for obsm key)")
    
    
    with open(os.path.join(out_dir, 'out', 'clustering_results.pkl'), 'rb') as f:
        clusters = pickle.load(f)

    # with open(os.path.join(out_dir, 'markers.pkl'), 'rb') as f:
    #     markers = pickle.load(f)

    # perform final_merging
    print(f"====================")
    print(f"Performing final merging...")
    print(f"clusters before merging: {len(clusters)}")

    clusters_after_merging, markers_after_merging = final_merge(
        adata=adata, 
        cluster_assignments=clusters, 
        # marker_genes=markers, # required for PCA, but optional if using a pre-computed latent space
        n_samples_per_clust=20, 
        random_seed=2024, 
        n_jobs = 30, # modify this to the number of cores you want to use
        return_markers_df=False, # return the pair-wise DE results for each cluster pair. If False (default), only return a set of markers (top 20 of up and down regulated genes in each pair comparison)
        final_merge_kwargs=merge_kwargs)

    print(f"Finished final merging. Total number of clusters after merging: {len(clusters_after_merging)}")

    with open(os.path.join(out_dir, 'out', "clustering_results_after_merging.pkl"), 'wb') as f:
        pickle.dump(clusters_after_merging, f)
    
    # determine datatype for markers_after_merging and save
    if markers_after_merging is None:
        print("Skipped calculating markers. Did not save markers.")
    elif isinstance(markers_after_merging, pd.DataFrame):
        markers_after_merging.to_csv(os.path.join(out_dir, 'out', 'markers_after_merging.csv'))
        print('Finished writing the pair-wise DE results to a .csv file.')
    else:
        with open(os.path.join(out_dir, 'out', 'markers_after_merging.pkl'), 'wb') as f:
            pickle.dump(markers_after_merging, f)
        print('Finished writing top markers to a .pkl file')
    
    # convert the clustering results to a .csv file
    n_cells = sum(len(i) for i in clusters_after_merging)
    cl = ['unknown']*n_cells
    for i in range(len(clusters_after_merging)):
        for j in clusters_after_merging[i]:
            cl[j] = i+1
    res = pd.DataFrame({'cl': cl}, index=adata.obs_names)
    res.to_csv(os.path.join(out_dir, 'out', 'clustering_results_after_merging.csv'))
    
    print(f"Finished writing clustering results to {out_dir}/out/")

if __name__ == "__main__":
    adata_path = sys.argv[1]
    latent_path = sys.argv[2]
    out_path = sys.argv[3]
    final_merge_kwargs = sys.argv[4]
    main(adata_path, latent_path, out_path, final_merge_kwargs)