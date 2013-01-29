#! /bin/sh

if [ $# -ne 0 ]; then
        echo "Param "$@" "
	ibrun -n $1 -o 0 $2 -d $3 $4

fi;
