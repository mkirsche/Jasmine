# Automated pipeline for alignment and SV calling in long-read datasets

### Installation:
Ensure that `snakemake` is installed.
Clone the repository via `git clone https://github.com/aganezov/EnsemblePLInternal.git` into a general location (e.g., `/path/to/pipelines`) on a computing cluster.

### Experiment run:
Determine and `cd` into experiment dedicated folder (e.g., `experiment`).
Create symlinks to all `snakefile` and `py` files in the `EnsemblePIInteral` master folder:
```bash
ln -s /path/to/pipelines/EnsemblePIInteral/*snakefile .
ln -s /path/to/pipelines/EnsemblePIInteral/*py .
``` 
Copy configuration `yaml` files:
```bash
cp /path/to/pipelines/EnsemblePIInteral/*yaml .
```
Set up targeted dataset inside `data.yaml` file.

Run snakemake (dry run) to see automated pipeline:
```bash
snakemake -s pipeline.snakefile -npr
``` 
If the dry run produces satisfactory results, run `snakemake` with production settings, including, possible cluster setup, multithreading, etc.

#### Example SLURM run:
```bash
snakemake -s pipeline.snakefile --latency-wait 200 -pr -j 20 --rerun-incomplete --cluster "sbatch --account={cluster.account} --partition={cluster.partition} --job-name={cluster.name} --nodes={cluster.nodes} --cpus-per-task={cluster.nCPUs} --time={cluster.time} --out={cluster.out} --err={cluster.err} --mem={cluster.mem_mb}M"
```
which ensures that
 * no more than 20 jobs (`-j`) are submitted at a time
 * any incomplete results from possible previous failed runs are regenerated (`--rerun-incomplete`)
 * sets up a SLURM submission setup, which in turn:
    * requests a single node of `parallel` partition per job,
    * with a 3 day time limit,
    * with 24G of RAM per node,
    * and with 
