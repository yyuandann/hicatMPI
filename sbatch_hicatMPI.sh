#!/bin/bash
#SBATCH --job-name=hicatMPI
#SBATCH --partition=celltypes
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=100:00:00
#SBATCH --mem=500G
#SBATCH --mail-user=dan.yuan@alleninstitute.org
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# Things to edit in this .sh file: 
# 1. Number of nodes and size to use. for over 1 million cells, I recommend using 11 nodes with each 500G memory. It will take 2 hours to finish. 
#  - 5 200G nodes are enough for 20k cells
#  - If HPC is busy, reduce the number of nodes. If requested n nodes, change the number of nodes for the workers to use in the last line of code to (n-1)
# 2. the absolute path to the h5ad file, scvi latent space (.csv file), and a output path where the results will be saved. The .h5ad file should have raw counts or normalized counts in adata.X. The algorithm will automatically detect if the counts are raw or normalized. If it is raw counts, it will do the normalization. This will be noted in the log file
# 3. the absolute path to the manager script and the worker script
# 4. the clustering parameters. The default parameters are set to match the bigcat default parameters. You can change them as needed.

scripts_dir="$1"
adata_path="$2"
latent_path="$3"
out_dir="$4"
clust_kwargs="$5"

manager_script="${scripts_dir}/iterative_clustering_mpi_manager.py"
worker_script="${scripts_dir}/iterative_clustering_mpi_worker.py"

source /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/etc/profile.d/conda.sh
conda activate /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/envs/tc

if [ ! -d "$out_dir" ]; then
    # If the directory doesn't exist, create it
    mkdir -p "$out_dir"
    echo "Directory '$out_dir' created."
fi
cd "$out_dir" # navigate to the output directory, where the log and out files will be saved

n_workers=$(( SLURM_NTASKS - 1 ))
echo "SLURM_NNODES=${SLURM_NNODES} SLURM_NTASKS=${SLURM_NTASKS} => manager=1 worker=${n_workers}"

# export PYTHONPATH=/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering:$PYTHONPATH
export PYTHONPATH="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering${PYTHONPATH:+:$PYTHONPATH}" # with "set -euo" the previous line would fail if PYTHONPATH was not set which is the case

time mpiexec \
    -n 1 sh -c "python \"$manager_script\" \"$adata_path\" \"$latent_path\" \"$out_dir\" \"$clust_kwargs\" > manager_output.log 2> manager_error.log" \
    : \
    -n "$n_workers" sh -c "python \"$worker_script\" \"$out_dir\" 2> worker_error.log" 