if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

if [ ! -d $BINDIR/Iris/src ]
then
    git submodule update --init --remote Iris
fi

javac -cp $BINDIR/Iris/src $BINDIR/src/*.java 

