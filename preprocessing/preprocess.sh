# Script for performing preprocessing for SV merging
# Converts duplications to insertion calls and uses Iris for refining insertion calls
# Assumes the external_scripts folder within Iris already has working scripts,
# and that samtools is on the user's path

BINDIR=`dirname $(readlink -f "$0")`

if [ $# -ne 3 ] && [ $# -ne 4 ]
then
  echo "USAGE: preprocess.sh file_list.txt [bam_file_list.txt] genome.fa out_file.vcf"
  exit
fi

REFINE=0

if [ $# -eq 3 ]
then
  VCFFILE=$1
  GENOME=$2
  OUTFILE=$3
  REFINE=0
fi

if [ $# -eq 4 ]
then
  VCFFILE=$1
  BAMFILE=$2
  GENOME=$3
  OUTFILE=$4
  REFINE=1
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

javac $BINDIR/*.java

VCFFILE_CONVERTED_DUPS=$VCFFILE'_duptoins.vcf'
java -cp $BINDIR DuplicationsToInsertions $VCFFILE $GENOME $VCFFILE_CONVERTED_DUPS

VCFFILE_SPEC=$VCFFILE_CONVERTED_DUPS'_spec'.vcf
java -cp $BINDIR MarkSpecificCalls $VCFFILE_CONVERTED_DUPS $VCFFILE_SPEC 10 30

if [ $REFINE -eq 1 ]
then
  $BINDIR/RefineInsertions/build.sh
  java -cp $BINDIR/RefineInsertions/ Iris genome_in=$GENOME vcf_in=$VCFFILE_SPEC reads_in=$BAMFILE vcf_out=$OUTFILE
else
  cp $VCFFILE_SPEC $OUTFILE
fi


