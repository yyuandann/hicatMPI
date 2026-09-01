from mpi4py import MPI
import queue
import pickle
import time
import copy
import os
import importlib
import shutil
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import re
import ast
import psutil
import anndata as ad
import scipy.sparse as sp
import warnings, logging
warnings.filterwarnings("ignore")          # silence pandas/numpy warning flood from the aligned merge
logging.disable(logging.WARNING)           # silence the clustering package's INFO/DEBUG flood (keeps ERROR/CRITICAL)
pd.options.mode.chained_assignment = None

# Which package provides the clustering implementation. Defaults to transcriptomic_clustering;
# set HICAT_PKG=pybigcat (with HICAT_PKG_PATH) to run the pipeline against the older package instead.
# The two sbatch scripts export both variables, so a run is self-describing.
_PKG = os.environ.get("HICAT_PKG", "transcriptomic_clustering")
_ic = importlib.import_module(f"{_PKG}.iterative_clustering")
onestep_clust, OnestepKwargs = _ic.onestep_clust, _ic.OnestepKwargs

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

MANAGER_RANK = 0

def preprocess_dict(dict_str):
# Remove comments using regex
    dict_str_cleaned = re.sub(r"#.*", "", dict_str)
    return dict_str_cleaned

def append_list_to_pkl(filepath, new_list):
    with open(filepath, "ab") as f:
        pickle.dump(new_list, f)

def append_markers(filepath, new_markers):
    with open(filepath, "rb") as f:
        markers = pickle.load(f)
        markers = markers | new_markers
    with open(filepath, "wb") as f:
        pickle.dump(markers, f)

def load_pkl(filepath):
    """Load all lists from a .pkl file into a list of lists."""
    results = []
    with open(filepath, "rb") as f: 
        while True:
            try:
                # Load the next object (a list) from the file
                cluster = pickle.load(f)
                results.append(cluster)
            except EOFError:
                # End of file reached, stop reading
                break
    return results

def get_max(X):

    if sp.issparse(X):
        max_val = X.data.max() if X.nnz else 0.0
        # sparse matrices are implicitly zero elsewhere
        max_val = max(max_val, 0.0)
    else:
        max_val = np.max(X)

    return max_val

def manager_job_queue(adata_path, latent_path, out_path, clust_kwargs,
                      random_seed=2024): # data is a tuple of anndata and key argumnents
    start = time.perf_counter()
        
    clust_kwargs = preprocess_dict(clust_kwargs)
    clust_kwargs = ast.literal_eval(clust_kwargs)

    # PCA-path blocks (means_vars/highly_variable/pca/filter_pcs/filter_known_modes) are used ONLY when
    # latent_component is None; on the latent/embedding path they are ignored. Use .get(...) so they can
    # be omitted entirely from clust_kwargs (default {} -> OnestepKwargs defaults are used).
    means_vars_kwargs = clust_kwargs.get('means_vars_kwargs', {})
    highly_variable_kwargs = clust_kwargs.get('highly_variable_kwargs', {})
    pca_kwargs = clust_kwargs.get('pca_kwargs', {})
    filter_pcs_kwargs = clust_kwargs.get('filter_pcs_kwargs', {})
    filter_known_modes_kwargs = clust_kwargs.get('filter_known_modes_kwargs', {})
    latent_kwargs = clust_kwargs.get('latent_kwargs', {})
    # split_size: recursion stop -- min cluster size to keep sub-clustering. Set inside clust_kwargs.
    split_size = clust_kwargs.get('split_size', 4)
    # save_markers: whether to save the union of split-time DE genes as markers_before_final_merge.csv.
    # Off by default: markers are properly computed on the FINAL clusters as a separate step
    # (de_all_pairs + marker queries); the split-time union is only needed as a diagnostic or as
    # input for a space='markers' final merge (runs without an integrated latent).
    save_markers = clust_kwargs.get('save_markers', False)
    cluster_louvain_kwargs = clust_kwargs.get('cluster_louvain_kwargs', {})
    merge_clusters_kwargs = clust_kwargs.get('merge_clusters_kwargs', {})

    # Build two parameter sets that differ ONLY in the merge DE score threshold:
    #   - clust_kwargs_top:    used by the manager's first (top-level) clustering
    #   - clust_kwargs_refine: used by the workers' recursive sub-clustering
    # merge_clusters_kwargs['thresholds']['score_thresh'] may be either a SCALAR (same threshold for
    # both levels) or a 2-element list [top, recursive] (mirrors R de.score.th = 300 top / 150 refine).
    merge_top = copy.deepcopy(merge_clusters_kwargs)
    merge_refine = copy.deepcopy(merge_clusters_kwargs)
    st = merge_clusters_kwargs.get('thresholds', {}).get('score_thresh')
    if isinstance(st, (list, tuple)):
        if len(st) != 2:
            raise ValueError(f"score_thresh must be a scalar or [top, recursive]; got {st}")
        top_th, rec_th = st[0], st[1]
    else:
        top_th = rec_th = st
    # downstream merge_clusters_by_de expects a scalar score_thresh, so resolve the list here
    merge_top['thresholds']['score_thresh'] = top_th
    merge_refine['thresholds']['score_thresh'] = rec_th

    # both share the same latent_kwargs dict, so setting latent_component below applies to both
    def _make_onestep_kwargs(merge_kw):
        return OnestepKwargs(
            means_vars_kwargs = means_vars_kwargs,
            highly_variable_kwargs = highly_variable_kwargs,
            pca_kwargs = pca_kwargs,
            filter_pcs_kwargs = filter_pcs_kwargs,
            filter_known_modes_kwargs = filter_known_modes_kwargs,
            latent_kwargs = latent_kwargs,
            cluster_louvain_kwargs = cluster_louvain_kwargs,
            merge_clusters_kwargs = merge_kw
        )

    clust_kwargs_top = _make_onestep_kwargs(merge_top)
    clust_kwargs_refine = _make_onestep_kwargs(merge_refine)

    # State the merge mode explicitly. Batch-aware merging fails SILENTLY if the kwarg is dropped or
    # misspelled -- the run completes normally, just pooled -- so a log line is the only thing that
    # distinguishes "ran batch-aware" from "quietly did not".
    _bam = merge_clusters_kwargs.get('batch_aware_merging')
    if _bam:
        merge_mode_msg = (f"BATCH-AWARE on obs[{_bam.get('batch_obs')!r}] "
                          f"(lfc_conservation_th={_bam.get('lfc_conservation_th', 0.7)}, "
                          f"per-batch overrides: {sorted(_bam.get('thresholds', {}) or {})})")
    else:
        merge_mode_msg = "POOLED (batch_aware_merging not set)"

    print(
        f"split_size={split_size} | "
        f"top-level score_thresh={merge_top['thresholds']['score_thresh']} | "
        f"recursive score_thresh={merge_refine['thresholds']['score_thresh']} | "
        f"merge={merge_mode_msg}",
        flush=True,
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
        print(f"Raw count data provided", flush=True)
        print(f"Normlazing total counts to 1e6...", flush=True)
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        new_max = get_max(adata.X)
        print(f"Finished normalization. max:{new_max}", flush=True)
    else:
        print(f"Normalized data provided. Skipped normalization", flush=True)
    
    print(f"Memory usage with adata in memory: {psutil.virtual_memory().percent}%", flush=True)

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

        adata.obsm['latent'] = np.asarray(latent)
        latent_kwargs['latent_component'] = 'latent'  # shared dict -> applies to both top & refine kwargs
        print(f"Finished reading in latent space: {latent.shape}", flush=True)

    else:
        print(f"latent_path is invalid or not provided. Using latent_kwargs['latent_component'] for clustering (None for PCA and a str for obsm key)", flush=True)
    
    # random_seed is now a function parameter (default 2024), set via CLI arg [5]
    print(f"random_seed={random_seed}", flush=True)

    # Batch-aware config vs the FULL data. This is the one place the complete modality set is
    # known, so typos are caught here; deeper in the recursion a branch may hold only some
    # modalities and merge_clusters treats overrides for the missing ones as inert.
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

    # Start each run from a clean slate. Leftover chunks from an interrupted run would either be
    # orphaned or mismatched to this run's cell partition; a missing tmp dir crashes the first
    # sub-cluster write. Wipe + recreate so a cancel-and-resubmit is always safe.
    tmp_dir_h5ads = os.path.join(out_path, 'tmp_h5ads')
    tmp_dir_idx   = os.path.join(out_path, 'tmp_idx')
    for _d in (tmp_dir_h5ads, tmp_dir_idx):
        shutil.rmtree(_d, ignore_errors=True)
        os.makedirs(_d, exist_ok=True)

    out_dir = os.path.join(out_path, "out")
    os.makedirs(out_dir, exist_ok=True)
    # clustering_results_tmp.pkl is written in append mode; drop any leftover so results aren't doubled.
    _leftover = os.path.join(out_dir, 'clustering_results_tmp.pkl')
    if os.path.exists(_leftover):
        os.remove(_leftover)

    # let the manager do the first clustering (top-level score threshold)
    clusters, markers = onestep_clust(adata, clust_kwargs_top, random_seed)
    sizes = [len(cluster) for cluster in clusters]
    print(f"Manager finished clustering {adata.shape[0]} cells into {len(clusters)} clusters of sizes {sizes}", flush=True)

    # save the markers from the 1st clustering (accumulator; only if save_markers)
    if save_markers:
        with open(os.path.join(out_dir, 'markers_before_final_merge.pkl'), 'wb') as f:
            pickle.dump(markers, f)
    
    # Initialize job queue, and put all subclusters into the queue
    adata_tmp_idx = 1
    job_queue = queue.Queue()
    if(len(clusters) == 1):
        append_list_to_pkl(os.path.join(out_dir, 'clustering_results_tmp.pkl'), clusters[0])
    else:
        for i in range(len(clusters)):
            if len(clusters[i]) < split_size:
                append_list_to_pkl(os.path.join(out_dir, 'clustering_results_tmp.pkl'), clusters[i])    
            else:
                idx = clusters[i]
                new_adata = adata[idx]
                new_adata_path = os.path.join(tmp_dir_h5ads, str(adata_tmp_idx)+'.h5ad')
                new_adata.write(new_adata_path)

                new_idx_path = os.path.join(tmp_dir_idx, str(adata_tmp_idx)+'.pkl')
                with open(new_idx_path, 'wb') as f:
                    pickle.dump(idx, f)

                new_task = (onestep_clust, (new_adata_path, clust_kwargs_refine, random_seed, new_idx_path))
                job_queue.put(new_task)

                adata_tmp_idx += 1

    # send the clusters from the first one-step clustering for further clustering 
    active_workers = 0

    for worker_rank in range(1, size):
        if not job_queue.empty():
            task = job_queue.get()
            comm.send(task, dest=worker_rank, tag=1)
            active_workers += 1
        else:
            comm.send(None, dest=worker_rank, tag=0) # terminate nodes that have no tasks
    print(f"Initiate {active_workers} activate workers", flush=True)

    # manage completed tasks from workers and assign new tasks
    while not job_queue.empty() or active_workers > 0:
        status = MPI.Status() # it can be defined outside of the while loop, does not matter
        clusters, new_markers, n_cells = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status) # receive any signals from any workers
        worker_rank = status.Get_source() # determine which worker just sent a message to the manager node 
        # tag = status.Get_tag() # tag==1: actual results, tag==0: acknowledgement of termination signal

        if save_markers:
            append_markers(os.path.join(out_dir, 'markers_before_final_merge.pkl'), new_markers) # internal accumulator; converted to .csv and removed at the end of the run
        p = psutil.virtual_memory().percent
        if(p > 90):
            print(f"Caution! Memory usage exceeds 90%: {p}%", flush=True)

        sizes = [len(cluster) for cluster in clusters]
        print(f"Worker {worker_rank} finished clustering {n_cells} cells into {len(clusters)} clusters of sizes {sizes}", flush=True)

        # save the finished clustering results and add new tasks to the queue
        if(len(clusters) == 1):
            append_list_to_pkl(os.path.join(out_dir, 'clustering_results_tmp.pkl'), clusters[0])

        else:
            for i in range(len(clusters)):
                if len(clusters[i]) < split_size:
                    append_list_to_pkl(os.path.join(out_dir, 'clustering_results_tmp.pkl'), clusters[i])

                else:
                    idx = clusters[i]
                    new_adata = adata[idx]
                    new_adata_path = os.path.join(tmp_dir_h5ads, str(adata_tmp_idx)+'.h5ad')
                    new_adata.write(new_adata_path)

                    new_idx_path = os.path.join(tmp_dir_idx, str(adata_tmp_idx)+'.pkl')
                    with open(new_idx_path, 'wb') as f:
                        pickle.dump(idx, f)

                    new_task = (onestep_clust, (new_adata_path, clust_kwargs_refine, random_seed, new_idx_path))
                    job_queue.put(new_task)

                    adata_tmp_idx += 1

        if not job_queue.empty():
            new_task = job_queue.get()
            comm.send(new_task, dest=worker_rank, tag=1)
        else:
            active_workers -= 1
            comm.send(None, dest=worker_rank, tag=0)
            print(f"Manger node sent None to worker node {worker_rank}. Currently {active_workers} active workers", flush=True)

    end = time.perf_counter()
    print(f"Finished all clustering tasks in {end - start:0.4f} seconds", flush=True)
    
    results = load_pkl(os.path.join(out_dir, 'clustering_results_tmp.pkl'))
    print(f"Total number of clusters: {len(results)}", flush=True)

    # convert the clustering results to a .csv file -- keyed by cell name, this is THE clustering
    # artifact downstream steps consume (index-based .pkl outputs were dropped: decoupled from the
    # adata they index into, they fail silently on a reordered/subset adata)
    n_cells = sum(len(i) for i in results)
    cl = ['unknown']*n_cells
    for i in range(len(results)):
        for j in results[i]:
            cl[j] = i+1
    res = pd.DataFrame({'cl': cl}, index=adata.obs_names)
    res.to_csv(os.path.join(out_dir, 'clusters_before_final_merge.csv'))

    # markers accumulated across all clustering jobs: write the user-facing copy as a plain
    # one-gene-per-line csv, then drop the pkl accumulator (internal state only)
    if save_markers:
        with open(os.path.join(out_dir, 'markers_before_final_merge.pkl'), 'rb') as f:
            _markers = pickle.load(f)
        pd.Series(sorted(_markers), name='gene').to_csv(
            os.path.join(out_dir, 'markers_before_final_merge.csv'), index=False)
        os.remove(os.path.join(out_dir, 'markers_before_final_merge.pkl'))
    else:
        print("save_markers=False: split-time markers not saved (compute markers on the final "
              "clusters with de_all_pairs as a separate step)", flush=True)

    # remove the tmp folders if empty:
    if os.path.exists(tmp_dir_h5ads) and len(os.listdir(tmp_dir_h5ads)) == 0:
        os.rmdir(tmp_dir_h5ads)
    else:
        print(f"tmp folder is not empty, sth went wrong", flush=True)
    
    if os.path.exists(tmp_dir_idx) and len(os.listdir(tmp_dir_idx)) == 0:
        os.rmdir(tmp_dir_idx)
    
    os.remove(os.path.join(out_dir, 'clustering_results_tmp.pkl'))

    print(f"Finished writing clustering results to {out_dir}", flush=True)
    print(f"Next step: run a final merge of clusters", flush=True)

if __name__ == "__main__":
    if rank == MANAGER_RANK:
        adata_path = sys.argv[1]
        latent_path = sys.argv[2]
        out_path = sys.argv[3]
        clust_kwargs = sys.argv[4]

        # Optional args (positional). Missing or the literal "None" -> use defaults:
        #   [5] random_seed : RNG seed (default 2024)
        # split_size and the top/recursive merge DE score thresholds are set INSIDE clust_kwargs
        # ('split_size', and merge_clusters_kwargs['thresholds']['score_thresh'] = [top, recursive]).
        def _opt(i):
            return sys.argv[i] if len(sys.argv) > i and sys.argv[i] not in ('', 'None') else None

        random_seed = int(_opt(5)) if _opt(5) is not None else 2024

        manager_job_queue(adata_path, latent_path, out_path, clust_kwargs, random_seed)