#!/bin/sh

host=linux-gfortran
prefix=/usr/local
debugsym=false
openmpsym=true

for arg in "$@"; do
    case "$arg" in
    --prefix=*)
        prefix=`echo $arg | sed 's/--prefix=//'`
        ;;

    --host=*)
        host=`echo $arg | sed 's/--host=//'`
        ;;

    --enable-debug)
        debugsym=true;;
    --disable-debug)
        debugsym=false;;

    --enable-openmp)
        openmpsym=true;;
    --disable-openmp)
        openmpsym=false;;

    --help)
        echo 'usage: ./configure [options]'
        echo 'options:'
        echo '  --prefix=<path>: installation prefix'
        echo '  --host=<system-name>: installation prefix'
	echo '  --list-systems: print valid system names'
        echo '  --enable-debug: include debug symbols'
        echo '  --disable-debug: do not include debug symbols'
	echo '  --enable-openmp: include openmp symbols'
	echo '  --disable-openmp: do not include openmp symbols'
        echo 'all invalid options are silently ignored'
        exit 0
        ;;
    esac
done

echo 'generating makefile ...'
echo "PREFIX = $prefix" >Makefile
echo "HOST = $host" >>Makefile
if $debugsym; then
    echo 'dbg = yes' >>Makefile
fi
if $openmpsym; then
    echo 'openmp = yes' >>Makefile
fi
cat Makefile.in >>Makefile
echo 'configuration complete, type make to build.'
