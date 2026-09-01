#!/bin/bash
#SBATCH --job-name=final_merging
#SBATCH --partition=celltypes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --time=100:00:00
# --mem is supplied by submit_pipeline_*.sh on the sbatch command line
# (command-line options override #SBATCH directives). --nodes stays 1: the final merge is single-node.

hicatMPI_scripts_dir="$1"
adata_path="$2"
latent_path="$3"
out_dir="$4"
final_merge_kwargs="$5"

wrapper_script="${hicatMPI_scripts_dir}/final_merging.py"

# Conda env: settable from submit_pipeline_*.sh via HICAT_CONDA_ENV (passed on the sbatch
# --export line); the default below is the fallback for direct invocation. Activation itself
# must happen HERE, in the job, so the environment is self-contained on the compute node.
HICAT_CONDA_ENV="${HICAT_CONDA_ENV:-/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/envs/tc}"
source /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/etc/profile.d/conda.sh
conda activate "$HICAT_CONDA_ENV"

cd "$out_dir"

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

time python "$wrapper_script" "$adata_path" "$latent_path" "$out_dir" "$final_merge_kwargs" > final_merge.log 2>&1
