#!/bin/bash
#SBATCH --job-name=hicatMPI
#SBATCH --partition=celltypes
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=100:00:00
# --nodes and --mem are supplied by submit_pipeline_*.sh on the sbatch command line
# (command-line options override #SBATCH directives).

set -euo pipefail

# Shared machinery -- do NOT edit per run. All user settings (paths, nodes, memory, conda env,
# package checkout, clustering/merge parameters) live in submit_pipeline_pooled.sh /
# submit_pipeline_batchaware.sh, which invoke this script and pass everything in.
# The worker count is derived from the allocation automatically (n_workers = SLURM_NTASKS - 1).

hicatMPI_scripts_dir="$1"
adata_path="$2"
latent_path="$3"
out_dir="$4"
clust_kwargs="$5"
# split_size (recursion stop) and the top/recursive DE score thresholds are set INSIDE clust_kwargs
# ('split_size', and merge_clusters_kwargs['thresholds']['score_thresh'] = [top, recursive]).
random_seed="${6:-None}"   # optional: RNG seed (default 2024 in manager); used for self-consistency runs

manager_script="${hicatMPI_scripts_dir}/iterative_clustering_mpi_manager.py"
worker_script="${hicatMPI_scripts_dir}/iterative_clustering_mpi_worker.py"

# Conda env: settable from submit_pipeline_*.sh via HICAT_CONDA_ENV (passed on the sbatch
# --export line); the default below is the fallback for direct invocation. Activation itself
# must happen HERE, in the job, so the environment is self-contained on the compute node.
HICAT_CONDA_ENV="${HICAT_CONDA_ENV:-/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/envs/tc}"
source /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/etc/profile.d/conda.sh
conda activate "$HICAT_CONDA_ENV"

if [ ! -d "$out_dir" ]; then
    # If the directory doesn't exist, create it
    mkdir -p "$out_dir"
    echo "Directory '$out_dir' created."
fi
cd "$out_dir" # navigate to the output directory, where the log and out files will be saved

n_workers=$(( SLURM_NTASKS - 1 ))
echo "SLURM_NNODES=${SLURM_NNODES} SLURM_NTASKS=${SLURM_NTASKS} => manager=1 worker=${n_workers}"

# Which clustering package to run against. Defaults to transcriptomic_clustering; export
# HICAT_PKG=pybigcat and HICAT_PKG_PATH=<...>/tool/pybigcat to use the older package instead.
export HICAT_PKG="${HICAT_PKG:-transcriptomic_clustering}"
export HICAT_PKG_PATH="${HICAT_PKG_PATH:-/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering}"
export PYTHONPATH="${HICAT_PKG_PATH}${PYTHONPATH:+:$PYTHONPATH}"   # ":+" so "set -u" tolerates an unset PYTHONPATH
echo "[hicatMPI] package=${HICAT_PKG} path=${HICAT_PKG_PATH}"
# Compute nodes vary in how old /lib64/libstdc++ is; some are older than this env needs
# (igraph -> libicuuc -> GLIBCXX_3.4.30), which makes the import fail on those nodes only.
# Prefer the env's own libraries so the job runs the same wherever SLURM places it.
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

time mpiexec \
    -n 1 sh -c "python \"$manager_script\" \"$adata_path\" \"$latent_path\" \"$out_dir\" \"$clust_kwargs\" \"$random_seed\" > manager_output.log 2> manager_error.log" \
    : \
    -n "$n_workers" sh -c "python \"$worker_script\" \"$out_dir\" 2> worker_error.log" 