#!/usr/bin/env bash
set -euo pipefail

# =====================================================================================
# hicatMPI pipeline — BATCH-AWARE merging (a pair merges only if every batch/modality votes to merge)
#
# Copy this file into your project directory, edit the paths and parameters, and run it.
# The copy then documents that run's exact configuration — keep it with the outputs.
#
# This is the multimodal counterpart of submit_pipeline_pooled.sh. The CLUSTERING is identical
# in both pipelines (embedding-based Louvain/Leiden -- R's harmonizing knn_joint clustering is
# deliberately not ported); what changes, and must change TOGETHER, is:
#   split_size 100          children keep enough cells PER PLATFORM (>= min.cells) to judge
#   batch_aware_merging     R merge_cl_multiple at every level and in the final merge
#
# The parameter values below are the WMB2 production settings validated against R
# (test_clustering/reports/10_harmonize_multimodal.md: identical merge partitions given
# identical inputs). Threshold values marked "dataset choice" should be revisited per dataset.
#
# INPUT CONTRACT
#   adata_path   .h5ad; adata.X may be RAW counts or ln(CPM+1)-normalized. Raw counts are
#                auto-detected (max value > 100) and normalized by the manager and final-merge
#                scripts using scanpy (sc.pp.normalize_total to 1e6 + sc.pp.log1p = natural log).
#                If you pre-normalize, it MUST be natural log, NOT log2 -- log2 data passes
#                the auto-detection silently but shifts every threshold below
#                (ln(2) = 0.6931472 corresponds to R's 1 in log2).
#                Batch-aware mode requires an IN-MEMORY AnnData -- size memory accordingly.
#                adata.obs must carry the batch column named by 'batch_obs'.
#   gene axis    prefer the UNION of the platforms' library gene lists: it keeps
#                single-platform evidence, which R counts (an INTERSECTION silently discards
#                it). A union axis REQUIRES 'genes_by_batch' (below) = each platform's library
#                list -- whoever built the union h5ad has these -- because the matrix itself
#                cannot distinguish 'never measured' from 'measured, zero counts'. An
#                intersection axis needs no such parameter (already exactly R).
#   latent_path  csv of the integrated embedding (e.g. scVI), first column = cell names
#                covering EVERY cell in the h5ad; overrides latent_component -> 'latent'
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
# Batch-aware merging REQUIRES the bigcat-aligned checkout.
tc_pkg_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering"

# Step 1 (clustering) always runs; Step 2 (global final merge) is optional.
run_final_merging=true              # false -> clustering only

clust_kwargs="{
    'split_size': 100,                           # R multimodal i_harmonize (harmonize.R:683)
    'save_markers': False,                       # markers are computed on the FINAL clusters as a
                                                 # separate step (de_all_pairs; see the tc examples).
                                                 # True saves the split-time DE-gene union as
                                                 # markers_before_final_merge.csv (diagnostic, or
                                                 # input for a space='markers' final merge)
    'latent_kwargs': {'latent_component': None}, # overridden to 'latent' by latent_path above
    'cluster_louvain_kwargs': {
        'k': 15,                                 # R knn_joint k
        'nn_measure': 'euclidean',               # R self.knn.method = Annoy.Euclidean
        'knn_method': 'pynndescent',             # library default; 'annoy' matches scrattch.bigcat
        'knn_seed': 1,                           # fixed -> seed-invariant graph (both backends)
        'weighting_method': 'jaccard_snn',       # SNN graph (R jaccard_big)
        'louvain_method': 'vtraag',              # Leiden, as R method='leiden'
        'resolution': 1.0,
        'jaccard_prune': 0.05,                   # same as the pooled pipeline (both use the same
                                                 # clustering); applied only to graphs > 50,000 cells
        'n_jobs': 30
    },
    'merge_clusters_kwargs': {
        'thresholds': {
            'q1_thresh': 0.4,                    # dataset choice -- WMB2 relaxed de.param
            'q2_thresh': None,                   #   (run_iter_cluster.R:54-70); R default is 0.5
            'qdiff_thresh': 0.7,
            'cluster_size_thresh': 10,           # R min.cells
            'padj_thresh': 0.01,
            'lfc_thresh': 0.6931472,             # ln(2) = R lfc.th=1 (log2)
            'low_thresh': 0.6931472,             # ln(2) = R low.th=1 (log2)
            'score_thresh': 100,                 # dataset choice -- R de.score.th, WMB2 relaxed
            'min_genes': 5,
            'merge_mode': 'aligned'
        },
        'k': 4,                                  # R compare.k
        'de_method': 'ebayes',                   # REQUIRED in batch-aware mode
        'batch_aware_merging': {                 # R merge_cl_multiple at EVERY recursion level --
            'batch_obs': 'platform',             # which is what knn_joint does (harmonize.R:640)
            'lfc_conservation_th': 0.7,          # R lfc.conservation.th
            'thresholds': {                      # dataset choice -- per-batch overrides (the full
                '10X_cells_v3':  {'q1_thresh': 0.4},   # WMB2 production setup, run_iter_cluster.R:54-70).
                '10X_cells_v2':  {'q1_thresh': 0.4},   # Unlisted batches inherit the shared thresholds
                '10X_nuclei_v3': {'q1_thresh': 0.3}    # above (so the two 0.4 lines are redundant with
            },                                   # q1_thresh 0.4 there -- kept to show the full setup).
                                                 # Nuclei sample less RNA, so their detection bar is
                                                 # lower. Any threshold can be overridden per batch.
            # 'conservation_lfc_th': 0.6931472,  # optional; defaults to thresholds['lfc_thresh'].
            #                                    # R hardcodes 1 (log2) = ln(2), not de.param\$lfc.th
            # 'genes_by_batch': {...},           # union gene axis only: {batch: [genes that batch's
            #                                    # library measures]} so each platform is judged on its
            #                                    # own gene list, as R does. Omit for an intersection
            #                                    # axis (every platform measures every gene there).
        }
    }
}"

# Final-merge params (used only if run_final_merging=true). With a latent_path the merge runs on
# the embedding and only these blocks are read. PCA path (no embedding) additionally needs:
#    'pca_kwargs': {'n_comps': 50, 'svd_solver': 'randomized'},
#    'filter_pcs_kwargs': {'known_components': None, 'similarity_threshold': 0.7, 'method': 'zscore', 'zth': 2, 'max_pcs': None},
#    'filter_known_modes_kwargs': {'known_modes': 'log2ngene', 'similarity_threshold': 0.7},
#    'project_kwargs': {},
final_merge_kwargs="{
    'merge_clusters_kwargs': {
        'thresholds': {
            'q1_thresh': 0.4,                    # dataset choice -- WMB2 relaxed de.param
            'q2_thresh': None,
            'qdiff_thresh': 0.7,
            'cluster_size_thresh': 10,           # R min.cells
            'padj_thresh': 0.01,
            'lfc_thresh': 0.6931472,             # ln(2) = R lfc.th=1 (log2)
            'low_thresh': 0.6931472,             # ln(2) = R low.th=1 (log2)
            'score_thresh': 100,                 # dataset choice -- R de.score.th, WMB2 relaxed
            'min_genes': 5,
            'merge_mode': 'aligned'
        },
        'k': 4,
        'de_method': 'ebayes',
        'n_markers': None,                       # None skips marker calc (faster)
        'batch_aware_merging': {                 # same keys and meaning as in clust_kwargs; makes
            'batch_obs': 'platform',             # the FINAL merge batch-aware too. The reduced space
            'lfc_conservation_th': 0.7,          # that shortlists pairs is unchanged (the embedding).
            'thresholds': {                      # the full WMB2 per-batch setup, as in clust_kwargs
                '10X_cells_v3':  {'q1_thresh': 0.4},
                '10X_cells_v2':  {'q1_thresh': 0.4},
                '10X_nuclei_v3': {'q1_thresh': 0.3}
            }
            # 'genes_by_batch': {...},           # as above, for a union gene axis
        }
    },
    'latent_kwargs': {'latent_component': None}
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

# Run:  chmod +x submit_pipeline_batchaware.sh && ./submit_pipeline_batchaware.sh
