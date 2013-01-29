
#!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    NP=16
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
#make

#-run the case       
echo Running with NP = $NP       

current_dir=$PWD
echo $current_dir

#data_dir=${current_dir/"/share/home/01523"/"/scratch/01523"}
data_dir=`echo $PWD|sed "s/share\/home/scratch/g"`


#llist=(.05)
llist=(1 .8 .7 .5 .3 .1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4)
#.5 .1 .2 .3)
#llist=(.05 
#llist=(.5 .3 .2 .004) 
#llist = (8e-4 5e-4)
#llist=(.3 5e-4)
#tstep=(1e-3 1e-3 1e-3 1e-3 1e-3 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2) 
#1e-3)
#1e-2 1e-2 1e-2 1e-2 1e-2)
key=3D_lin
queue=normal
llist=(.1 5e-2 1e-2 5e-3 5e-4)
#key=fast_g64
#rm status_$key.log

allruns=()
for lval in ${llist[@]}
do
    current_dir=$data_dir/data_${key}_${lval}
    echo $current_dir
  

    
    cp q3_simp.cxx $current_dir/physics_code.cxx.ref #universal convension for now

done
