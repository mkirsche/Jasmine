# Thriver

Tool for merging structural variants.

It's still in development, so while it works, it doesn't yet have the full set of features.

## Instructions for building

./build.sh

## Instructions for running

java -cp src Main

Running this without parameters will provide a menu describing the necessary parameters.

## Optional Input INFO Fields

Some parameters are custimizable on a per-variant basis.  In these cases, the value for each variant can be specified by assing an extra INFO field to each VCF entry containing the value.  A list of these is below:

* THRIVER_DIST - the maximum distance other variants can be to merge with this one
* THRIVER_ID - the minimum sequence identity other insertions need to have to this variant to be able to merge with it


## Methods Description

The document describing the method can be found here: https://docs.google.com/document/d/1UNYLisKn8OpPRdYGoxQMN174baGTaxsG5CuuR1QuCV0/edit?usp=sharing
