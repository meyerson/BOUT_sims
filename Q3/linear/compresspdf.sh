#!/bin/bash

if [ $# -ne 0 ]; then
        echo "Param "$@" "
	gs -dBATCH -dNOPAUSE -dSAFER -sDEVICE=jpeg -dJPEGQ=5 -r200x200 -sOutputFile=tempjpg-%03d.jpeg $1 
	convert tempjpg*.jpeg small_$1
	rm tempjpg*

fi