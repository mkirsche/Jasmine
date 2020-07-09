if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

git submodule update $BINDIR/Iris

irisjar=''
if [ $# -eq 1 ]
then
    echo "Iris jar: " $1
    irisjar=$1
else
    $BINDIR/Iris/build_jar.sh
    irisjar=$BINDIR/Iris/iris.jar
fi

cd $BINDIR/src
javac -cp $irisjar *.java
jar -c -e Main -f jasmine.jar *.class
mv jasmine.jar $BINDIR
cd $BINDIR
