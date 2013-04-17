 #!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    NP=8
    MPIEXEC="mpirun -np"
fi  
# Usage: scriptname -options
# Note: das(-) necessary

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



make

rm status.log
#rm run.log



sim_key='3D'
codekey='SOLblob'
inpkey='3D'

current_dir=$PWD
data_dir=/tmp/${codekey}

current_dir=$data_dir/data_${sim_key}
echo $current_dir

#rm -r $current_dir
mkdir -p $current_dir

rm -r $PWD/data
  

cp ${codekey}.cxx   $current_dir/physics_code.cxx.ref
cp BOUT_${inpkey}.inp $current_dir/BOUT.inp
 
echo $MPIEXEC
    
ln -s $current_dir $PWD/data
echo $current_dir >> status.log
    #echo $current_dir


$MPIEXEC $NP ./${codekey} -d $current_dir 
  