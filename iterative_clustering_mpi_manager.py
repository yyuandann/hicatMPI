from mpi4py import MPI
import queue
import pickle
import time
import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import re
import ast
import psutil
import anndata as ad

from transcriptomic_clustering.iterative_clustering import onestep_clust, OnestepKwargs

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

def manager_job_queue(adata_path, latent_path, out_path, clust_kwargs): # data is a tuple of anndata and key argumnents
    start = time.perf_counter()
        
    clust_kwargs = preprocess_dict(clust_kwargs)
    clust_kwargs = ast.literal_eval(clust_kwargs)

    means_vars_kwargs = clust_kwargs['means_vars_kwargs']
    highly_variable_kwargs = clust_kwargs['highly_variable_kwargs']
    pca_kwargs = clust_kwargs['pca_kwargs']
    filter_pcs_kwargs = clust_kwargs['filter_pcs_kwargs']
    filter_known_modes_kwargs = clust_kwargs['filter_known_modes_kwargs']
    latent_kwargs = clust_kwargs['latent_kwargs']
    cluster_louvain_kwargs = clust_kwargs['cluster_louvain_kwargs']
    merge_clusters_kwargs = clust_kwargs['merge_clusters_kwargs']

    clust_kwargs = OnestepKwargs(
        means_vars_kwargs = means_vars_kwargs,
        highly_variable_kwargs = highly_variable_kwargs,
        pca_kwargs = pca_kwargs,
        filter_pcs_kwargs = filter_pcs_kwargs,
        filter_known_modes_kwargs = filter_known_modes_kwargs,
        latent_kwargs = latent_kwargs,
        cluster_louvain_kwargs = cluster_louvain_kwargs,
        merge_clusters_kwargs = merge_clusters_kwargs
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
        print(f"Raw count data provided", flush=True)
        print(f"Normlazing total counts to 1e6...", flush=True)
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        print(f"Finished normalization. max:{np.max(adata.X)}", flush=True)
    else:
        print(f"Normalized data provided. Skipped normalization", flush=True)
    
    print(f"Memory usage with adata in memory: {psutil.virtual_memory().percent}%", flush=True)

    # determine if latent_path is valid
    if os.path.exists(latent_path):
        latent = pd.read_csv(latent_path, index_col=0)
        latent = latent.loc[adata.obs_names]
        adata.obsm['latent'] = np.asarray(latent)
        clust_kwargs.latent_kwargs['latent_component'] = 'latent'
        print(f"Finished reading in latent space: {latent.shape}", flush=True)
    else:
        print(f"latent_path is invalid or not provided. Using latent_kwargs['latent_component'] for clustering (None for PCA and a str for obsm key)", flush=True)
    
    min_samples = 4
    random_seed = 2024

    tmp_dir_h5ads =  os.path.join(out_path,'tmp_h5ads')
    if not os.path.exists(tmp_dir_h5ads):
        os.makedirs(tmp_dir_h5ads)

    tmp_dir_idx =  os.path.join(out_path,'tmp_idx')
    if not os.path.exists(tmp_dir_idx):
        os.makedirs(tmp_dir_idx)

    out_dir = os.path.join(out_path, "out")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # let the manager do the first clustering
    clusters, markers = onestep_clust(adata, clust_kwargs, random_seed)
    sizes = [len(cluster) for cluster in clusters]
    print(f"Manager finished clustering {adata.shape[0]} cells into {len(clusters)} clusters of sizes {sizes}", flush=True)

    # save the markers from the 1st clustering
    with open(os.path.join(out_dir, 'markers.pkl'), 'wb') as f:
        pickle.dump(markers, f)
    
    # Initialize job queue, and put all subclusters into the queue
    adata_tmp_idx = 1
    job_queue = queue.Queue()
    if(len(clusters) == 1):
        append_list_to_pkl(os.path.join(out_dir, 'clustering_results_tmp.pkl'), clusters[0])
    else:
        for i in range(len(clusters)):
            if len(clusters[i]) < min_samples:
                append_list_to_pkl(os.path.join(out_dir, 'clustering_results_tmp.pkl'), clusters[i])    
            else:
                idx = clusters[i]
                new_adata = adata[idx]
                new_adata_path = os.path.join(tmp_dir_h5ads, str(adata_tmp_idx)+'.h5ad')
                new_adata.write(new_adata_path)

                new_idx_path = os.path.join(tmp_dir_idx, str(adata_tmp_idx)+'.pkl')
                with open(new_idx_path, 'wb') as f:
                    pickle.dump(idx, f)

                new_task = (onestep_clust, (new_adata_path, clust_kwargs, random_seed, new_idx_path))
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

        append_markers(os.path.join(out_dir, 'markers.pkl'), new_markers) # writes to the markers.pkl file
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
                if len(clusters[i]) < min_samples:
                    append_list_to_pkl(os.path.join(out_dir, 'clustering_results_tmp.pkl'), clusters[i])

                else:
                    idx = clusters[i]
                    new_adata = adata[idx]
                    new_adata_path = os.path.join(tmp_dir_h5ads, str(adata_tmp_idx)+'.h5ad')
                    new_adata.write(new_adata_path)

                    new_idx_path = os.path.join(tmp_dir_idx, str(adata_tmp_idx)+'.pkl')
                    with open(new_idx_path, 'wb') as f:
                        pickle.dump(idx, f)

                    new_task = (onestep_clust, (new_adata_path, clust_kwargs, random_seed, new_idx_path))
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

    with open(os.path.join(out_dir, "clustering_results.pkl"), 'wb') as f:
        pickle.dump(results, f)
    print(f"Total number of clusters: {len(results)}", flush=True)

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
        manager_job_queue(adata_path, latent_path, out_path, clust_kwargs)