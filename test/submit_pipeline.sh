#!/usr/bin/env bash
set -euo pipefail

scripts_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/hicatMPI_v2"
adata_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/QCed_gex_byNeigh_final_cleaned/marmoset_TH-EPI-Glut.zarr"
latent_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/projectedLatent/marmoset_TH-EPI-Glut.csv"
out_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/hicatMPI_v2/test"

n_nodes_mpi=5
mem=200G

# hard-coded kwargs here
# the following are the default parameters (matching the defaults in scratch.bigcat) for the clustering. change them as needed.
clust_kwargs="{
    'means_vars_kwargs': {
        'low_thresh': 0.6931472, # lowest value required for a gene to pass filtering. set to 1 originally, change to 0.6931472 to match to the bigcat default
        'min_cells': 4 # minimum number of cells expressed required for a gene to pass filtering
    },
    'highly_variable_kwargs': {
        'max_genes': 4000 # originally 3000, change to 4000 to match to the bigcat default
    },
    'pca_kwargs': {
        'cell_select': 30000, # originally 500000 cells
        'n_comps': 50,
        'svd_solver': 'randomized'
    },
    'filter_pcs_kwargs': {
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore', # or elbow
        'zth': 2,
        'max_pcs': None,
    },
    # if not using known_modes, set filter_known_modes_kwargs to None or an empty dict. Only applies for PCA
    'filter_known_modes_kwargs': {
        'known_modes': 'log2ngene', 
        'similarity_threshold': 0.7
    },
    ## !!NEW!! Original method: "PCA", allows the user to select any obsm latent space such as "X_scVI" for leiden clustering.
    'latent_kwargs': {
        'latent_component': None # None (default to run PCA) or obsm key such as "X_scVI"
    },
    'cluster_louvain_kwargs': {
        'k': 15, # number of nn, originally 150, change to 15 to match to the bigcat default
        'nn_measure': 'euclidean',
        'knn_method': 'annoy',
        'louvain_method': 'taynaud', #'vtraag',
        'weighting_method': 'jaccard',
        'n_jobs': 30, # cpus, originally 8
        'resolution': 1.0 # resolution of louvain for taynaud method
    },
    'merge_clusters_kwargs': {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, ## originally uses 50, change to 10 to match to the bigcat default
            'qdiff_thresh': 0.7, 
            'padj_thresh': 0.05, 
            'lfc_thresh': 1, # log2 fold change threshold for DE genes
            'score_thresh': 100, # originally uses 200, change to 100 to match to the bigcat default
            'low_thresh': 0.6931472, # originally uses 1 # applied to log2(cpm+1) to determine if a gene is expressed or not, change to 0.6931472 to match to the bigcat default
            'min_genes': 5
        },
        'k': 4, # number of nn for de merge, originaly 2, change to 4 to match to the bigcat default
        'de_method': 'ebayes'
    }
}"

final_merge_kwargs="{
    'pca_kwargs':{
        # 'cell_select': 30000, # should not use set this for final merging, as we need to sample from each cluster if computing PCA
        'n_comps': 50,
        'svd_solver': 'randomized'
    },
    'filter_pcs_kwargs': {
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore', 
        'zth': 2,
        'max_pcs': None
    },
    'filter_known_modes_kwargs': {
        'known_modes': 'log2ngene', 
        'similarity_threshold': 0.7
    },
    'project_kwargs': {},
    'merge_clusters_kwargs': {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, 
            'qdiff_thresh': 0.7, 
            'padj_thresh': 0.05, 
            'lfc_thresh': 1, 
            'score_thresh': 100, 
            'low_thresh': 0.6931472, 
            'min_genes': 5
        },
        'k': 4,
        'de_method': 'ebayes',
        'n_markers': None, # if set to None, it will bypass the marker calculation step, which speeds up things
    },
    'latent_kwargs': { 
        'latent_component': None # None or a obsm in adata. if None: default is running pca, else use the latent_component in adata.obsm
    }
}"

cd $out_dir

# submit the hicatMPI job
mpi_submit_out=$(sbatch --nodes="$n_nodes_mpi" --mem="$mem" "${scripts_dir}/sbatch_hicatMPI.sh" \
  "$scripts_dir" "$adata_path" "$latent_path" "$out_dir" "$clust_kwargs")

mpi_jobid=$(echo "$mpi_submit_out" | awk '{print $4}')
echo "Submitted mpi.sh as job $mpi_jobid with nodes=$n_nodes_mpi"

# submit merge job after MPI success
merge_submit_out=$(sbatch --mem="$mem" --dependency=afterok:"$mpi_jobid" "${scripts_dir}/sbatch_finalMerging.sh" \
  "$scripts_dir" "$adata_path" "$latent_path" "$out_dir" "$final_merge_kwargs")

merge_jobid=$(echo "$merge_submit_out" | awk '{print $4}')
echo "Submitted final_merging.sh as job $merge_jobid (will run after $mpi_jobid)"

# submit the full job in bash:
# chmod +x submit_pipeline.sh
# ./submit_pipeline.sh
