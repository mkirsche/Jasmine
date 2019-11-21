# Script for performing preprocessing for SV merging
# Reverts back insertions which were originally duplications

BINDIR=`dirname $(readlink -f "$0")`

if [ $# -ne 2 ]
then
  echo "USAGE: postprocess.sh svs.vcf out.vcf"
  exit
fi
VCFFILE=$1
OUTFILE=$2

if [ ! -r $VCFFILE ]
then
  echo "Cannot read sv vcf file: $VCFFILE"
  exit
fi

javac $BINDIR/*.java

java -cp $BINDIR InsertionsToDuplications $VCFFILE $OUTFILE
