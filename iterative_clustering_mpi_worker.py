from mpi4py import MPI
import anndata as ad
import os
import pickle
import traceback
import sys
import warnings, logging
warnings.filterwarnings("ignore")          # silence pandas/numpy warning flood from the aligned merge
logging.disable(logging.WARNING)           # silence the clustering package's INFO/DEBUG flood (keeps ERROR/CRITICAL)
import pandas as _pd; _pd.options.mode.chained_assignment = None

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

MANAGER_RANK = 0

def worker_job(out_path):
    # Create a directory for worker logs if it doesn't exist
    os.makedirs(os.path.join(out_path, "worker_logs"), exist_ok = True)

    with open(f"{out_path}/worker_logs/worker_{rank}_output.log", 'w') as out_file:
        #  open(f"{out_path}/worker_logs/worker_{rank}_error.log", 'w') as err_file: # did not work
        # Redirect stdout to output file
        sys.stdout = out_file
        # sys.stderr = err_file # did not work, will just write to the slurm output file

        while True: # this is necessary to keep the worker running, otherwise it will just process just one task before exit
            # Receive a job from the manager
            task = comm.recv(source=MANAGER_RANK, tag=MPI.ANY_TAG)
        
            if task is None:
                print(f"Worker {rank} received termination signal.", flush=True)
                break  

            try:
                func, args = task
                adata_path, kwargs, random_seed, idx_path = args

                with open(idx_path, 'rb') as f:
                    idx = pickle.load(f)
                
                # total number of idices in idx (a list of lists)
                idx_size = len(idx)

                adata = ad.read_h5ad(adata_path)
                n = adata.shape[0]

                if idx_size != n:
                    raise ValueError(f"Mismatch in number of indices. idx: {idx_size},  adata: {n}")

                clusters, markers = func(adata, kwargs, random_seed)

                # convert the indices back to the original indices
                original_clusters = []
                for cluster in clusters:
                    original_clusters.append([idx[i] for i in cluster])

                # Signal completion to the manager node (need to save to disk and pass paths instead of objects as well?)
                comm.send((original_clusters, markers, n), dest=MANAGER_RANK) 

                os.remove(adata_path)
                os.remove(idx_path)
            
            except Exception as e:
                error_msg = f"Worker {rank} error: {str(e)}\n{traceback.format_exc()}"
                # err_file.write(error_msg + "\n")
                # err_file.flush()
                print(f"Error occurred: {error_msg}", flush=True)
                raise
            

if __name__ == "__main__":
    out_path = sys.argv[1]
    worker_job(out_path)
