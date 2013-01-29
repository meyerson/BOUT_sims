#!/bin/bash

module add mencoder

llist=(1e-1 5e-3)
key=3d
i=0

data_dir=`echo $PWD|sed "s/share\/home/scratch/g"`
for lval in ${llist[@]}
do
    current_dir=$data_dir/data_${key}_${lval}
    qsub -q serial -pe 1way 16 -N "movie"_$lval  -cwd -V -m e -A BOUT++_startup -l h_rt=4:00:00 -j y makemovie.py $current_dir $lval

    echo "$((i++))"
done