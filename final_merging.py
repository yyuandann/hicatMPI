import os
import importlib
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import re
import ast
import time
import scipy.sparse as sp

# Which package provides the clustering implementation. Defaults to transcriptomic_clustering;
# set HICAT_PKG=pybigcat (with HICAT_PKG_PATH) to run the pipeline against the older package instead.
# The two sbatch scripts export both variables, so a run is self-describing.
_PKG = os.environ.get("HICAT_PKG", "transcriptomic_clustering")
_fm = importlib.import_module(f"{_PKG}.final_merging")
final_merge, FinalMergeKwargs = _fm.final_merge, _fm.FinalMergeKwargs

def get_max(X):
    if sp.issparse(X):
        max_val = X.data.max() if X.nnz else 0.0
        # sparse matrices are implicitly zero elsewhere
        max_val = max(max_val, 0.0)
    else:
        max_val = np.max(X)
    return max_val

def main(adata_path, latent_path, out_dir, final_merge_kwargs): # data is a tuple of anndata and key argumnents
    start = time.perf_counter()

    def preprocess_dict(dict_str):
    # Remove comments using regex
        dict_str_cleaned = re.sub(r"#.*", "", dict_str)
        return dict_str_cleaned
    
    final_merge_kwargs = preprocess_dict(final_merge_kwargs)
    final_merge_kwargs = ast.literal_eval(final_merge_kwargs)

    # The PCA-path blocks (pca / filter_pcs / filter_known_modes / project) are used ONLY when the
    # merge runs on space='pca' (no embedding). With a latent supplied they are ignored, so .get()
    # lets submit files omit them entirely ({} -> FinalMergeKwargs defaults).
    pca_kwargs = final_merge_kwargs.get('pca_kwargs', {})
    filter_pcs_kwargs = final_merge_kwargs.get('filter_pcs_kwargs', {})
    filter_known_modes_kwargs = final_merge_kwargs.get('filter_known_modes_kwargs', {})
    project_kwargs = final_merge_kwargs.get('project_kwargs', {})
    merge_clusters_kwargs = final_merge_kwargs['merge_clusters_kwargs']
    latent_kwargs = final_merge_kwargs.get('latent_kwargs', {})
    
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

    adata.obs_names = adata.obs_names.astype(str)

    if get_max(adata.X) > 100:
        print(f"Raw count data provided")
        print(f"Normlazing total counts to 1e6...")
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        new_max = get_max(adata.X)
        print(f"Finished normalization. max:{new_max}")
    else:
        print(f"Normalized data provided")

    # determine if latent_path is valid
    if os.path.exists(latent_path):
        latent = pd.read_csv(latent_path, index_col=0)
        latent.index = latent.index.astype(str)

        n_adata = adata.n_obs
        n_latent = latent.shape[0]
        n_shared = latent.index.isin(adata.obs_names).sum()

        print(
            f"Shared cells between adata and latent: {n_shared:,} | "
            f"removed from adata: {n_adata - n_shared:,} | "
            f"removed from latent: {n_latent - n_shared:,}",
            flush=True,
        )

        # keep only shared cells, preserve adata order
        keep = adata.obs_names[adata.obs_names.isin(latent.index)]
        adata = adata[keep].copy()
        latent = latent.loc[keep].copy()

        # latent = latent.loc[adata.obs_names]
        adata.obsm['latent'] = np.asarray(latent)
        merge_kwargs.latent_kwargs['latent_component'] = 'latent'
        print(f"Finished reading in latent space: {latent.shape}")
    else:
        print(f"scvi_path is invalid or not provided. Using latent_kwargs['latent_component'] for clustering (None for PCA and a str for obsm key)")
    
    
    # Load the clustering by CELL NAME (clusters_before_final_merge.csv), never by raw index:
    # a name-keyed join is order-independent and fails loudly if the adata and the clustering
    # have been decoupled (reordered, subset, or simply the wrong file).
    cl_df = pd.read_csv(os.path.join(out_dir, 'out', 'clusters_before_final_merge.csv'), index_col=0)
    cl_df.index = cl_df.index.astype(str)
    labels = cl_df['cl'].reindex(adata.obs_names)
    n_unlabeled = int(labels.isna().sum())
    if n_unlabeled:
        examples = list(adata.obs_names[labels.isna()][:3])
        raise ValueError(
            f"{n_unlabeled} of {adata.n_obs} cells in the adata have no label in "
            f"clusters_before_final_merge.csv (e.g. {examples}). The adata and the clustering "
            f"are out of sync -- refusing to merge on misaligned inputs.")
    extra = cl_df.index.difference(adata.obs_names)
    if len(extra):
        raise ValueError(
            f"{len(extra)} cells in clusters_before_final_merge.csv are absent from the adata "
            f"(e.g. {list(extra[:3])}). The adata and the clustering are out of sync -- "
            f"refusing to merge on misaligned inputs.")
    codes, uniques = pd.factorize(labels.values)
    # lists, not arrays: merge_clusters mutates assignments with list.extend()
    clusters = [np.flatnonzero(codes == i).tolist() for i in range(len(uniques))]

    # Batch-aware config vs the FULL data: catch modality-name typos here, where the complete
    # modality set is known, and report which modalities fall back to the shared thresholds.
    _bam = merge_clusters_kwargs.get('batch_aware_merging')
    if _bam:
        batch_obs = _bam.get('batch_obs')
        if batch_obs not in adata.obs:
            raise ValueError(f"batch_aware_merging['batch_obs']={batch_obs!r} is not a column of "
                             f"adata.obs (columns: {list(adata.obs.columns)})")
        modalities = set(map(str, pd.unique(adata.obs[batch_obs])))
        overrides = _bam.get('thresholds', {}) or {}
        unknown = sorted(set(map(str, overrides)) - modalities)
        if unknown:
            raise ValueError(
                f"batch_aware_merging['thresholds'] names modalities not present in the data: "
                f"{unknown}. Modalities in adata.obs[{batch_obs!r}]: {sorted(modalities)}")
        base_th = merge_clusters_kwargs.get('thresholds', {})
        for m in sorted(modalities):
            ov = overrides.get(m, {})
            if ov:
                print(f"[batch-aware] modality {m!r}: overrides {ov} (unlisted thresholds inherit "
                      f"the shared values)", flush=True)
            else:
                shown = {k: base_th.get(k) for k in ('q1_thresh', 'score_thresh') if k in base_th}
                print(f"[batch-aware] modality {m!r}: no per-modality override -- using the shared "
                      f"thresholds (e.g. {shown})", flush=True)

    # perform final_merging
    print(f"====================")
    print(f"Performing the final merge...")
    print(f"clusters before the final merge: {len(clusters)}")

    clusters_after_merging, markers_after_merging = final_merge(
        adata=adata, 
        cluster_assignments=clusters, 
        # marker_genes=markers, # required for the 'pca'/'markers' spaces, unused on the latent path.
        #                       # To supply: markers = set(pd.read_csv(os.path.join(out_dir, 'out',
        #                       #     'markers_before_final_merge.csv'))['gene'])
        n_samples_per_clust=20, 
        random_seed=2024, 
        n_jobs = 30, # modify this to the number of cores you want to use
        return_markers_df=False, # return the pair-wise DE results for each cluster pair. If False (default), only return a set of markers (top 20 of up and down regulated genes in each pair comparison)
        final_merge_kwargs=merge_kwargs)

    print(f"Finished the final merge. Total number of clusters after the final merge: {len(clusters_after_merging)}")

    # determine datatype for markers_after_merging and save
    if markers_after_merging is None:
        print("Skipped calculating markers. Did not save markers.")
    elif isinstance(markers_after_merging, pd.DataFrame):
        markers_after_merging.to_csv(os.path.join(out_dir, 'out', 'markers_after_final_merge.csv'))
        print('Finished writing the pair-wise DE results to a .csv file.')
    else:
        pd.Series(sorted(markers_after_merging), name='gene').to_csv(
            os.path.join(out_dir, 'out', 'markers_after_final_merge.csv'), index=False)
        print('Finished writing top markers to a .csv file (one gene per line)')
    
    # convert the clustering results to a .csv file
    n_cells = sum(len(i) for i in clusters_after_merging)
    cl = ['unknown']*n_cells
    for i in range(len(clusters_after_merging)):
        for j in clusters_after_merging[i]:
            cl[j] = i+1
    res = pd.DataFrame({'cl': cl}, index=adata.obs_names)
    res.to_csv(os.path.join(out_dir, 'out', 'clusters_after_final_merge.csv'))
    
    print(f"Finished writing clustering results to {out_dir}/out/")

if __name__ == "__main__":
    adata_path = sys.argv[1]
    latent_path = sys.argv[2]
    out_path = sys.argv[3]
    final_merge_kwargs = sys.argv[4]
    main(adata_path, latent_path, out_path, final_merge_kwargs)