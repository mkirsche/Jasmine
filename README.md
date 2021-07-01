[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/jasminesv/README.html)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=jasminesv)




# Jasmine

JASMINE: Jointly Accurate Sv Merging with Intersample Network Edges

Version 1.1.2

This tool is used to merge structural variants (SVs) across samples.  Each sample has a number of SV calls, consisting of position information (chromosome, start, end, length), type and strand information, and a number of other values.  Jasmine represents the set of all SVs across samples as a network, and uses a modified minimum spanning forest algorithm to determine the best way of merging the variants such that each merged variants represents a set of analogous variants occurring in different samples.


## Conda Installation

The recommended installation method is through [bioconda](https://bioconda.github.io/).

Conda Installation command (typically takes under a minute to install):

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda install jasminesv
```


## Instructions for building from source

When running Jasmine, one of the preprocessing options is to run Iris, a tool which refines the sequences and breakpoints of insertions in datasets with high-error reads.  Iris depends on samtools, minimap2, and racon by default, which can be installed separately and either added to your path or pointed to with the `iris_args` parameter.  Once these dependencies are installed (or if running Jasmine without Iris preprocessing), Jasmine can be built with the following command:

```
path_to_jasmine_repo/build_jar.sh
```


## Instructions for running

After building the jar file, Jasmine can be run with the executable file `jasmine`, which will be in the main folder of this repository if building from source, or in the condabin folder if installed through conda.  Running it with no parameters will print a usage menu describing the required and optional arguments.


## Demo Dataset
To run Jasmine on HiFi data from the HG002 trio, run the following commands (typically takes about a minute to download and under five minutes to run on a modern desktop):
```
wget https://bx.bio.jhu.edu/data/jasmine/HG002Trio/UnmergedVCFs/HG002vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf .
wget https://bx.bio.jhu.edu/data/jasmine/HG002Trio/UnmergedVCFs/HG003vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf .
wget https://bx.bio.jhu.edu/data/jasmine/HG002Trio/UnmergedVCFs/HG004vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf .
wget https://bx.bio.jhu.edu/data/jasmine/HG002Trio/HG002Trio_HiFi.merged.vcf .
ls *vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf > filelist.txt
jasmine file_list=filelist.txt out_file=merged.vcf
jasmine --dup_to_ins --postprocess_only out_file=merged.vcf
```

The output of merged.vcf should then exactly match the contents of HG002Trio_HiFi.merged.vcf.


## Optimized SV Inference Pipeline

Jasmine is offered as standalone software and will accurately merge SV calls from any SV callers, including short-read callers. However, if calling SVs from genomic long reads (PacBio CLR, PacBio HiFi, or Oxford Nanopore), for best results, we recommend using the following optimized pipeline to obtain population-scale SV calls from FASTQ files.  This pipeline is provided as a [Snakemake pipeline](https://github.com/mkirsche/Jasmine/tree/master/pipeline). 


![Jasmine SV Inference Pipeline](https://github.com/mkirsche/Jasmine/blob/master/pipeline/pipelineoverview.svg)



## IGV visualization module

Jasmine also includes a module for automating the creation of [IGV](http://software.broadinstitute.org/software/igv/) screenshots of variants of interest.  It can be run through the `igv_jasmine` executable file.  Running it with no parameters will print a usage menu describing the required and optional arguments, and it requires at minimum the following:
- BAM files from which variants were called in each sample
- The reference genome
- The merged VCF file, or a BED file with regions of interest

Running this module creates a folder which will store IGV screenshots for each variant (optionally filtered based on the command line parameters), and populates that folder with a .bat file, a script which can be run through IGV by selecting Tools -> Run Batch Script and navigating to the file.  After running this script, the folder containing the .bat file will also include images of the regions surrounding each variant of interest.


## User Manual

The user manual with detailed information about input/output files and command line arguments can be found here: https://github.com/mkirsche/Jasmine/wiki/Jasmine-User-Manual

