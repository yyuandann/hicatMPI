#!/bin/bash
#SBATCH --job-name=final_merging
#SBATCH --partition=celltypes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --time=100:00:00
#SBATCH --mem=500G
#SBATCH --mail-user=dan.yuan@alleninstitute.org
#SBATCH --mail-type=END,FAIL

scripts_dir="$1"
adata_path="$2"
latent_path="$3"
out_dir="$4"
final_merge_kwargs="$5"

wrapper_script="${scripts_dir}/final_merging.py"

source /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/etc/profile.d/conda.sh
conda activate /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/envs/tc

cd "$out_dir"

export PYTHONPATH=/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering:$PYTHONPATH

time python "$wrapper_script" "$adata_path" "$latent_path" "$out_dir" "$final_merge_kwargs" > final_merge.log 2>&1
