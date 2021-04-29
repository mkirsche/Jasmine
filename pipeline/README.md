# Automated pipeline for alignment and SV calling in long-read datasets

## Installation
Ensure that `snakemake` is installed.
Clone the repository via `git clone https://github.com/aganezov/EnsemblePLInternal.git` into a general location (e.g., `/path/to/pipelines`) on a computing cluster.

## Experiment run
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

### Example SLURM run
```bash
snakemake -s pipeline.snakefile --latency-wait 200 -pr -j 20 --rerun-incomplete --cluster "sbatch --account={cluster.account} --partition={cluster.partition} --job-name={cluster.name} --nodes={cluster.nodes} --cpus-per-task={cluster.nCPUs} --time={cluster.time} --out={cluster.out} --err={cluster.err} --mem={cluster.mem_mb}M"
```
which ensures that
 * no more than 20 jobs (`-j`) are submitted at a time
 * any incomplete results from possible previous failed runs are regenerated (`--rerun-incomplete`)
 * sets up a SLURM submission setup, which in turn requests a single node of `parallel` partition per job with:
    * a 3 day time limit,
    * 24G of RAM per node,

## Pipeline Overview

1. Align reads in each sample with [Winnowmap](https://github.com/marbl/Winnowmap) with the recommended parameters for the type of sequencing data being used
2. Call SVs in each sample with [sniffles](https://github.com/fritzsedlazeck/Sniffles) using sensitive parameters and report all supporting reads: `sniffles -m <alignments> -v <vcf> --threads <threads> --min_support 2 --max_distance 50  --min_length 20 --num_reads_report -1`
3. Convert duplications to insertions temporarily for breakpoint refinement and better cross-sample comparison: `jasmine --dup_to_ins --preprocess_only vcf_filelist=<vcf> --comma_filelist`
4. Refine SVs in each sample with [Iris](https://github.com/mkirsche/Iris/)
5. Normalize SV types in each sample: `jasmine --preprocess_only --pre_normalize --comma_filelist file_list=<vcf>`
6. Mark high-confidence callset (high-specificity callset) in each sample: `jasmine file_list=<vcf> --comma_filelist --preprocess_only --mark_specific spec_reads=<min(10, 25% average coverage)> spec_len=30`
7. Remove duplicate calls in each sample: `jasmine file_list=<vcf> max_dist=200 --allow_intrasample out_file=<outputvcf> --comma_filelist --nonlinear_dist`
8. Generate a list of all finalized per-sample VCF files (txt file, one per line)
9. Merge SVs across samples: `jasmine file_list=<vcflist> out_file=<outputmergedvcf>`
10. Convert insertions back to duplications: `jasmine --dup_to_ins --postprocess_only out_file=<mergedvcf>`
11. Remove low-confidence or imprecise calls: `cat <mergedvcf> | grep -v 'IMPRECISE;' | grep -v 'IS_SPECIFIC=0'`
