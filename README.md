# Jasmine

JASMINE: Jointly Accurate Sv Merging with Intersample Network Edges

Version 1.0.10

This tool is used to merge structural variants (SVs) across samples.  Each sample has a number of SV calls, consisting of position information (chromosome, start, end, length), type and strand information, and a number of other values.  Jasmine represents the set of all SVs across samples as a network, and uses a modified minimum spanning forest algorithm to determine the best way of merging the variants such that each merged variants represents a set of analogous variants occurring in different samples.


## Conda Installation

The recommended installation method is through [bioconda](https://bioconda.github.io/).

Conda Installation command:

```
conda install jasminesv
```


## Instructions for building from source

When running Jasmine, one of the preprocessing options is to run Iris, a tool which refines the sequences and breakpoints of insertions in datasets with high-error reads.  Iris depends on samtools, minimap2, and racon by default, which can be installed separately and either added to your path or pointed to with the `iris_args` parameter.  Once these dependencies are installed (or if running Jasmine without Iris preprocessing), Jasmine can be built with the following command:

```
path_to_jasmine_repo/build_jar.sh
```


## Instructions for running

After building the jar file, Jasmine can be run with the executable file `jasmine`, which will be in the main folder of this repository if building from source, or in the condabin folder if installed through conda.  Running it with no parameters will print a usage menu describing the required and optional arguments.


## IGV visualization module

Jasmine also includes a module for automating the creation of [IGV](http://software.broadinstitute.org/software/igv/) screenshots of variants of interest.  It can be run through the `igv_jasmine` executable file.  Running it with no parameters will print a usage menu describing the required and optional arguments, and it requires at minimum the following:
- BAM files from which variants were called in each sample
- The reference genome
- The merged VCF file, or a BED file with regions of interest

Running this module creates a folder which will store IGV screenshots for each variant (optionally filtered based on the command line parameters), and populates that folder with a .bat file, a script which can be run through IGV by selecting Tools -> Run Batch Script and navigating to the file.  After running this script, the folder containing the .bat file will also include images of the regions surrounding each variant of interest.

## User Manual

The user manual with detailed information about input/output files and command line arguments can be found here: https://github.com/mkirsche/Jasmine/wiki/Jasmine-User-Manual

