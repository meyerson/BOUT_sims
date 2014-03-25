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

zlist=(2)
data_dir='/tmp/hw'

for zval in ${zlist[@]}
do

  current_dir=$data_dir/data
  mkdir -p $current_dir
  # mkdir data_${zval}
  # ln -s data_${zval} data
  sed "s/Zeff = 128.0/Zeff = ${zval}.0/g" BOUT.inp > $current_dir/BOUT.inp
  if [ $zval -lt 128 ]
      then
      # reduce time-step. At large times these cases produce noise
      sed "s/TIMESTEP = 5e3/TIMESTEP = 1e3/g" $current_dir/BOUT.inp > $current_dir/tmp
      mv -f $current_dir/tmp $current_dir/BOUT.inp
  fi
      
  $MPIEXEC $NP ./hw -d $current_dir 
  rm -f data
done

#-check the result
#idl runidl.pro
