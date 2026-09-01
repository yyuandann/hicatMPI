#!/usr/bin/env bash
set -euo pipefail

# =====================================================================================
# hicatMPI pipeline — POOLED merging (modality-agnostic; the original pipeline)
#
# Copy this file into your project directory, edit the paths and parameters, and run it.
# The copy then documents that run's exact configuration — keep it with the outputs.
#
# To do batch-aware merging, use submit_pipeline_batchaware.sh instead: the merge
# strategy and split_size must change TOGETHER, so the two cases are separate
# self-contained files rather than comment toggles in one. (The clustering itself is
# identical in both pipelines. Pooled merging remains a valid choice for multi-modality
# data too -- it simply does not consult the modality when testing pairs.)
#
# INPUT CONTRACT
#   adata_path   .h5ad; adata.X may be RAW counts or ln(CPM+1)-normalized. Raw counts are
#                auto-detected (max value > 100) and normalized by the manager and final-merge
#                scripts using scanpy (sc.pp.normalize_total to 1e6 + sc.pp.log1p = natural log).
#                If you pre-normalize, it MUST be natural log, NOT log2 -- log2 data passes
#                the auto-detection silently but shifts every threshold below
#                (ln(2) = 0.6931472 corresponds to R's 1 in log2).
#   latent_path  csv of the embedding, first column = cell names covering EVERY cell in
#                the h5ad; when set it overrides latent_component -> 'latent'
#   out_dir      an 'out/' folder is created here
# =====================================================================================

hicatMPI_scripts_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/hicatMPI"
adata_path="/path/to/adata.h5ad"
latent_path="/path/to/latent.csv"
out_dir="/path/to/out_dir"
n_nodes_mpi=10                      # sizing guide: ~5 nodes x 200G for 20k cells;
mem=500G                            # 11 nodes x 500G cluster ~1M cells in ~2 h. One node is
                                    # the manager; the rest are workers.
mail_user=""                        # your email for job END/FAIL notifications; empty = no emails
conda_env="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/envs/tc"
                                    # conda env the jobs activate (must have the tc package deps)

# Which transcriptomic_clustering checkout the jobs import (passed explicitly to both sbatch
# scripts below; the import name is always 'transcriptomic_clustering', only the location varies).
tc_pkg_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering"

# Step 1 (clustering) always runs; Step 2 (global final merge) is optional.
run_final_merging=true              # false -> clustering only

# Clustering params. Defaults match R scrattch.bigcat's SINGLE-MODALITY path
# (iter_clust / iter_clust_big).
clust_kwargs="{
    'split_size': 10,                            # R single-modality iter_clust (cluster.R:367)
    'save_markers': False,                       # markers are computed on the FINAL clusters as a
                                                 # separate step (de_all_pairs; see the tc examples).
                                                 # True saves the split-time DE-gene union as
                                                 # markers_before_final_merge.csv (diagnostic, or
                                                 # input for a space='markers' final merge)
    'latent_kwargs': {'latent_component': None}, # None=PCA; a latent_path (above) overrides to the embedding
    'cluster_louvain_kwargs': {
        'k': 15,                                 # cell KNN neighbors
        'nn_measure': 'euclidean',
        'knn_method': 'pynndescent',             # library default; 'annoy' matches scrattch.bigcat --
                                                 # if switching to annoy, also set 'annoy_trees': 50
                                                 # (R's fixed value; the package default derives fewer)
        'knn_seed': 1,                           # fixed, decoupled from random_seed -> seed-invariant graph (both backends)
        'weighting_method': 'jaccard_snn',       # SNN graph (R jaccard_big)
        'louvain_method': 'vtraag',              # Leiden
        'resolution': 1.0,
        'jaccard_prune': 0.05,                   # R single-modality value (cluster_big.R:165);
                                                 # applied only to graphs > 50,000 cells, as R does
        'n_jobs': 30
    },
    'merge_clusters_kwargs': {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'qdiff_thresh': 0.7,
            'cluster_size_thresh': 10,           # R min.cells
            'padj_thresh': 0.01,
            'lfc_thresh': 0.6931472,             # ln(2) = R lfc.th=1 (log2)
            'low_thresh': 0.6931472,             # ln(2) = R low.th=1 (log2)
            'score_thresh': [300, 150],          # [top, recursive] (R de.score.th); scalar = both levels
            'min_genes': 5,
            'merge_mode': 'aligned'              # R-faithful merge ('fast' = quicker, coarser)
        },
        'k': 4,                                  # candidate merge neighbors (R merge_cl_big k)
        'de_method': 'ebayes'
    }
}"

# PCA path (no embedding): set latent_component=None and add these blocks into clust_kwargs:
#    'means_vars_kwargs': {'low_thresh': 0.6931472, 'min_cells': 4},
#    'highly_variable_kwargs': {'max_genes': 4000},
#    'pca_kwargs': {'cell_select': 30000, 'n_comps': 50, 'svd_solver': 'randomized'},
#    'filter_pcs_kwargs': {'known_components': None, 'similarity_threshold': 0.7, 'method': 'zscore', 'zth': 2, 'max_pcs': None},
#    'filter_known_modes_kwargs': {'known_modes': 'log2ngene', 'similarity_threshold': 0.7},

# Final-merge params (used only if run_final_merging=true). With a latent_path the merge runs on
# the embedding and only these blocks are read. PCA path (no embedding) additionally needs:
#    'pca_kwargs': {'n_comps': 50, 'svd_solver': 'randomized'},
#    'filter_pcs_kwargs': {'known_components': None, 'similarity_threshold': 0.7, 'method': 'zscore', 'zth': 2, 'max_pcs': None},
#    'filter_known_modes_kwargs': {'known_modes': 'log2ngene', 'similarity_threshold': 0.7},
#    'project_kwargs': {},
final_merge_kwargs="{
    'merge_clusters_kwargs': {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'qdiff_thresh': 0.7,
            'cluster_size_thresh': 10,           # R min.cells
            'padj_thresh': 0.01,
            'lfc_thresh': 0.6931472,             # ln(2) = R lfc.th=1 (log2)
            'low_thresh': 0.6931472,             # ln(2) = R low.th=1 (log2)
            'score_thresh': 100,                 # R de.score.th for the final merge
            'min_genes': 5,
            'merge_mode': 'aligned'
        },
        'k': 4,
        'de_method': 'ebayes',
        'n_markers': None                        # None skips marker calc (faster)
    },
    'latent_kwargs': {'latent_component': None}  # overridden to 'latent' by latent_path above
}"

# ---- submit: Step 1 clustering (always), Step 2 final merge (optional) ----
pkg_export="--export=ALL,HICAT_PKG_PATH=${tc_pkg_path},HICAT_CONDA_ENV=${conda_env}"
mail_opts=""
[ -n "$mail_user" ] && mail_opts="--mail-user=${mail_user} --mail-type=END,FAIL"
mpi_submit_out=$(sbatch $pkg_export $mail_opts --nodes="$n_nodes_mpi" --mem="$mem" "${hicatMPI_scripts_dir}/sbatch_hicatMPI.sh" \
  "$hicatMPI_scripts_dir" "$adata_path" "$latent_path" "$out_dir" "$clust_kwargs")
mpi_jobid=$(echo "$mpi_submit_out" | awk '{print $4}')
echo "Submitted clustering as job $mpi_jobid (nodes=$n_nodes_mpi, package=$tc_pkg_path)"

if [ "$run_final_merging" = true ]; then
  merge_submit_out=$(sbatch $pkg_export $mail_opts --mem="$mem" --dependency=afterok:"$mpi_jobid" "${hicatMPI_scripts_dir}/sbatch_finalMerging.sh" \
    "$hicatMPI_scripts_dir" "$adata_path" "$latent_path" "$out_dir" "$final_merge_kwargs")
  merge_jobid=$(echo "$merge_submit_out" | awk '{print $4}')
  echo "Submitted final_merging as job $merge_jobid (afterok:$mpi_jobid)"
else
  echo "run_final_merging=false -> clustering only."
fi

# Run:  chmod +x submit_pipeline_pooled.sh && ./submit_pipeline_pooled.sh
