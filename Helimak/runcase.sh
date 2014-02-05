#!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    NP=4
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


NOUTS=(300 20)
tstep=(2e0)
llist=(1e1)



i=0


key='hlmk_2D'


rm status_${key}.log

for lval in ${llist[@]}
do
  mkdir data_${lval}
  ln -s data_${lval} data
  
  current_dir=$data_dir/data_${key}_${lval}
  echo $current_dir
    
  #rm -r $current_dir
  mkdir -p $current_dir
  
  rm -r $PWD/data_${lval}

  cp ${key}.cxx   $current_dir/2fluid.cxx.ref
  cp ${key}.cxx   $current_dir/hlmk.cxx.ref
  cp ${key}.cxx   $current_dir/physics_code.cxx.ref

  cp ${key}.cxx   $PWD/data_${lval}/hlmk.cxx.ref
  #cp $data_dir/data_bz_${key}_${lval}/*restart* $current_dir

  sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT_${key}.inp > temp.inp
  #sed "s/ZMAX = 1/ZMAX = 1/g" BOUT.inp > temp.inp
  sed "s/NOUT = 100/NOUT = ${NOUTS[$i]}/g" temp.inp > temp2.inp
  #sed "s/NOUT = 100/NOUT = 100/g" temp.inp > temp2.inp
#sed "s/TIMESTEP = 5e2/TIMESTEP =  ${tstep[$i]}/g" temp2.inp > $current_dir/BOUT.inp
  sed "s/TIMESTEP = 5e2/TIMESTEP =  ${tstep[$i]}/g" temp2.inp > $current_dir/BOUT.inp
 
  echo "$((i++))"  
  
  $MPIEXEC $NP ./${key} -d $current_dir 
  ln -s $current_dir $PWD/data_${lval}
  echo $current_dir >> status_${key}.log
  #ibrun -n $NP -o 0  ./2fluid 
  #wait
  rm -f data

  python2.7 ./makemovie2.py $current_dir/ hlmk 1 ${NOUTS[$i]} ${NOUTS[$i]}-1
done
