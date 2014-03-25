#!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    NP=8
    MPIEXEC="mpirun -np"

fi  
# Usage: scriptname -options
# Note: dash (-) necessary

while getopts ":n:np" Option
do
  case $Option in
    n ) MPIEXEC="mpirun -np";NP=$OPTARG;;
    * ) ;;   # DEFAULT
  esac
done

#-compile/build local executable
make

#-run the case       
echo Running with NP = $NP       

rm -rf data*


current_dir=$PWD
data_dir='/tmp/hlmk'


NOUTS=(500)
tstep=(4e0)
llist=(5e0)



i=0


key='hlmk_2D'


rm status_${key}.log

for lval in ${llist[@]}
do
  
  
  current_dir=$data_dir/data_${key}_${lval}_te
  mkdir -p $current_dir
  # cp $PWD/restart/*  $current_dir
  # ls $current_dir
 
    
  #rm -r $current_dir
  

  cp ${key}.cxx   $current_dir/2fluid.cxx.ref
  cp ${key}.cxx   $current_dir/hlmk.cxx.ref
  cp ${key}.cxx   $current_dir/physics_code.cxx.ref

  sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT_${key}.inp > temp.inp
  #sed "s/ZMAX = 1/ZMAX = 1/g" BOUT.inp > temp.inp
  sed "s/NOUT = 100/NOUT = ${NOUTS[$i]}/g" temp.inp > temp2.inp
  sed "s/NXPE = 4/NXPE = ${NP}/g" temp2.inp > temp.inp
  sed "s/TIMESTEP = 5e2/TIMESTEP =  ${tstep[$i]}/g" temp.inp > $current_dir/BOUT.inp
 
  echo "$((i++))"  
  
  $MPIEXEC $NP ./${key} -d $current_dir 
  ln -s $current_dir $PWD/data_${lval}
  echo $current_dir >> status_${key}.log
  #ibrun -n $NP -o 0  ./2fluid 
  #wait
  rm -f data

  python2.7 ./makemovie2.py $current_dir/ hlmk 1 ${NOUTS[$i]} ${NOUTS[$i]}-1
done
