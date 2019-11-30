# Jasmine

JASMINE: Jointly Accurate Sv Merging with Intersample Network Edges

This tool is used to merge structural variants (SVs) across samples.  Each sample has a number of SV calls, consisting of positon information (chromosome, start, end, length), as well as type and strand information, and a number of other values.  Jasmine represents the set of all SVs across samples as a network, and uses a modified minimum spanning forest algorithm to determine the best way of merging the variants such that the merged variants represent analogous variants occurring in different samples.

## Instructions for building

./build.sh

## Instructions for running

java -cp PATHTOJASMINE/src:PATHTOJASMINE/Iris/src Main

Running this without any arguments will provide a menu describing the necessary and optional parameters.

## User Manual

The user manual with detailed information about input/output files and command line arguments can be found here: https://github.com/mkirsche/Jasmine/wiki/Jasmine-User-Manual

