# Script for performing preprocessing for SV merging
# Converts duplications to insertion calls and uses Iris for refining insertion calls
# Assumes the external_scripts folder within Iris already has working scripts,
# and that samtools is on the user's path

BINDIR=`dirname $(readlink -f "$0")`
VCFFILE=$1
BAMFILE=$2
GENOME=$3
OUTFILE=$4

if [ $# -ne 4 ]
then
  echo "USAGE: preprocess.sh svs.vcf longreads.bam genome.fa out.vcf"
  echo ""
  exit
fi

if [ ! -r $GENOME ]
then
  echo "Cannot read genome file: $GENOME"
  exit
fi

if [ ! -r $BAMFILE ]
then
  echo "Cannot read long reads bam file: $BAMFILE"
  exit
fi

if [ ! -r $VCFFILE ]
then
  echo "Cannot read sv vcf file: $VCFFILE"
  exit
fi

VCFFILE_CONVERTED_DUPS='duptoins_'$VCFFILE
javac $BINDIR/*.java
java -cp $BINDIR DuplicationsToInsertions $VCFFILE $GENOME $VCFFILE_CONVERTED_DUPS

VCFFILE_SPEC='spec_'$VCFFILE_CONVERTED_DUPS
java -cp $BINDIR MarkSpecificCalls $VCFFILE_CONVERTED_DUPS $VCFFILE_SPEC 10 30

$BINDIR/RefineInsertions/build.sh
java -cp $BINDIR/RefineInsertions/ Iris genome_in=$GENOME vcf_in=$VCFFILE_SPEC reads_in=$BAMFILE vcf_out=$OUTFILE



