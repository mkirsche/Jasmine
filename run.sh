BINDIR=`dirname $(readlink -f "$0")` # Where this bash script lives
WORKINGDIR=`pwd` # Where the script is run from

if [ $# -ne 3 ] && [ $# -ne 4 ]
then
  echo "USAGE: run.sh file_list.txt [bam_file_list.txt] genome.fa out_file.vcf"
  exit
fi

REFINE=0

if [ $# -eq 3 ]
then
  FILELIST=$1
  GENOME=$2
  OUTFILE=$3
  REFINE=0
fi

if [ $# -eq 4 ]
then
  FILELIST=$1
  BAMFILELIST=$2
  GENOME=$3
  OUTFILE=$4
  REFINE=1
fi

PREPROCESSED_FILELIST=$FILELIST'_preprocessed.txt'
if [ -r $PREPROCESSED_FILELIST ]
then
  rm $PREPROCESSED_FILELIST
fi

for i in `cat $FILELIST`
do
  echo 'Preprocessing ' $i
  PREPROCESSED_VCF=$i'_preprocessed'.vcf
  if [ $REFINE -eq 1 ]
  then
    $BINDIR/preprocessing/preprocess.sh $i $BAMFILE $GENOME $PREPROCESSED_VCF
  else
    $BINDIR/preprocessing/preprocess.sh $i $GENOME $PREPROCESSED_VCF
  fi

  echo $PREPROCESSED_VCF >> $PREPROCESSED_FILELIST
done

THRIVER_OUT=$WORKINGDIR/thriver_out.vcf

javac $BINDIR/src/*.java
java -cp $BINDIR/src Main file_list=$PREPROCESSED_FILELIST out_file=$THRIVER_OUT

$BINDIR/postprocessing/postprocess.sh $THRIVER_OUT $OUTFILE
