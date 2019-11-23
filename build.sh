if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

WORKINGDIR=`pwd`

cd $BINDIR
git submodule update --init --recursive
cd $WORKINGDIR
$BINDIR/Iris/build.sh

javac -cp $BINDIR/Iris/src $BINDIR/src/*.java
