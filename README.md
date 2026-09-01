# Accelerating *transcriptomic_clustering* with Distributed Computing


This method implements a **dynamic, asynchronous clustering algorithm** using the **Message-Passing Interface (MPI)** to distribute clustering tasks across multiple HPC nodes. It provides the same clustering results as the transcriptomic_clustering package (except that nearby cluster ids no longer imply similarity), but is optimized for distributed computating using MPI, significantly reducing runtime for large datasets. In v2.0.0, the final merging step is included in the same run.

## Quick Start

1. Pick the example that matches your data — **the entire user-facing configuration lives in this
   one file**:
   - **`submit_pipeline_pooled.sh`** — **pooled merging**: the original pipeline, ignoring
     batch/modality when merging clusters (one DE test over all cells; `split_size=10`).
   - **`submit_pipeline_batchaware.sh`** — **batch-aware merging**: every batch/modality votes,
     and a pair merges only if all of them vote to merge (R's `merge_cl_multiple`, at every level
     and in the final merge; `split_size=100` so children keep enough cells per batch to judge).
     The clustering itself is identical in both pipelines; the merge strategy and `split_size`
     belong together — that is why the two cases are separate files rather than comment toggles.
2. **Copy** the example into your project directory and edit it: input paths, thresholds, nodes and
   memory. The copy documents that run's exact configuration — keep it next to the outputs.
   - Input contract (both pipelines): the `.h5ad` must be **ln(CPM+1)** normalized (natural log) and
     the latent csv must cover every cell. Batch-aware mode additionally needs the AnnData to fit in
     memory, `adata.obs` to carry the batch column, and a union gene axis preferred (single-platform
     evidence is kept) with `genes_by_batch` = the platforms' library lists; an intersection axis
     needs no such parameter (see the file's header).
3. Make it executable and run it — this submits the clustering job and (optionally) the final-merge
   job with a dependency on it:
   ```bash
   chmod +x submit_pipeline_batchaware.sh
   ./submit_pipeline_batchaware.sh
   ```
4. Find the final outputs
   - The clustering results before and after the final merge are saved in clusters_before_final_merge.csv and clusters_after_final_merge.csv, respectively. For both .csv files, the index contains cell names and the cl column contains cluster IDs. All outputs are keyed by cell name (no index-based .pkl outputs): the final merge joins the clustering to the adata by cell name and refuses to run if they are out of sync.
   - Markers are computed on the final clusters as a **separate step** (`de_all_pairs` + marker queries; see the transcriptomic_clustering examples). Optionally, `'save_markers': True` in clust_kwargs saves the split-time DE-gene union as markers_before_final_merge.csv (a diagnostic, or input for a `space='markers'` final merge).
   - Refer to manager_output.log and final_merge.log for details on how clusters were split in each clustering job, and the number of clusters before and after final merging.

The batch-aware configuration ships with the WMB2 production values as a worked example; it was
validated against R's `merge_cl_multiple` (identical merge partitions given identical inputs) — see
[the alignment reports](https://github.com/yyuandann/transcriptomic_clustering/tree/bigcat-alignment/docs/tests_aligning_to_bigcat)
(report 10) for the validation and for which thresholds are dataset choices.

## A note on determinism

Which cells group together is deterministic — but the cluster IDs do not represent cluster
similarity or mean anything (unlike R's, where nearby IDs imply similar clusters); they only
reflect the order in which the clustering finished.

## Background
The transcriptomic clustering Python package that uses a **scVI latent space** can be found here: [transcriptomic_clustering](https://github.com/AllenInstitute/transcriptomic_clustering/tree/dev) (the `dev` branch; the former `hmba/tc_latent` branch was merged into it via PRs #131/#144 and deleted). It is the Python version of the R package [scrattch.hicat](https://github.com/AllenInstitute/scrattch.hicat), both of which perform clustering recursively (depth-first search). The recursive approach can take significant time for large datasets. For instance, clustering 1 million cells can take ~2 days.

## Distributed Clustering Model
This method replaces the depth-first search (DFS) recursive method with a **dynamic, asynchronous manager-worker model**:
- **Manager Node**: Oversees the clustering process, distributing jobs to worker nodes and managing the queue of tasks.
- **Worker Nodes**: Independently perform clustering tasks and, once a job finishes, immediately taking a job from the queue and performing clustering again.
- **Asynchronous Task Distribution**: The system doesn’t wait for all clustering jobs at the same hierarchy to finish before moving on. Instead, as soon as a worker node completes its current job, it directly moves to the next pending task in the queue.

### How It Works

![Process Illustration](images/mpiTC.jpeg)

- **Step 1**: The manager node performs the initial clustering task. Once that finishes, the manager node appends all subsequent clustering jobs for each cluster into the job queue. For each available worker node, the manager assigns a job from the queue.
- **Step 2**: Each worker node performs its assigned clustering task.
- **Step 3**: Once a worker finishes, it sends the result to the manager and immediately takes a pending job from the queue. The manager evaluate the clusters:
  - for clusters that cannot be further clustered (clustered into 1 cluster), the results are added to the final results.
  - for clusters that can be furtehr clustered (clustered into > 1 clusters), the subclusters are added to the queue for further clustering.
- **Step 4**: The process continues until the job queue is empty and all nodes have terminated.

## Benchmarking Results
### Run time comparison
The table below shows the comparison of run time between the transcriptomic_clustering package and hicatMPI.
| Number of Cells | transcriptomic_clustering | hicatMPI (This Method)|
|-----------------|---------------------------|-----------------------|
| 10,000          | 12 minutes                | 9 minutes             |
| 80,000          | 63 minutes                | 18 minutes            |
| 1,000,000       | 47 hours 41 minutes       | 2 hours 15 minutes    |

- Clustering of 3.3 million cells using hicatMPI was completed on 10 nodes (500 GB each) in 7 hours.

### Clustering results comparison
The confusion matrices below shows the clustering results from the transcriptomic_clustering package (y axis) and hicatMPI (x axis) for clustering 10k cells (left), 78k cells (middle), and 1 million cells (right). The 1-1 correspondence confirms that the same clustering results are obtained from both the original transcriptomic_clustering package and this MPI-based implementation, given the same data input and clustering parameters.

<img src="images/10k.png" width="250"/> <img src="images/78k.png" width="237"/> <img src="images/1m.png" width="240"/>
