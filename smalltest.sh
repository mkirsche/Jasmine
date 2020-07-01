if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

$BINDIR/build_no_iris.sh
$BINDIR/run.sh file_list=$BINDIR/test_data/a.vcf,$BINDIR/test_data/b.vcf out_file=$BINDIR/test_data/merged.vcf --comma_filelist
 
myout=$BINDIR/test_data/merged.vcf
correctout=$BINDIR/test_data/c.vcf

diff -w $myout $correctout >/dev/null;REPLY=$?
echo ''
if [ ${REPLY} -eq 0 ]
then
  echo '### TEST SUCCEEDED ###'
else
  echo '### TEST FAILED ###'
  diff -w $myout $correctout
fi
