
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
#llist=(.1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4)
#llist=(1e-3 5e-4 5e-2 1e-2 5e-3)
#llist=(.1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4)
llist=(1 .8 .7 .5 .3 .1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4)
# .1 .2 .3)
tstep=(5e-3 5e-3 5e-3)
#1e-3 1e-3 1e-3)
key=re_eig64
#rm status_$key.log

make

allruns=()
for lval in ${llist[@]}
do
    current_dir=$data_dir/data_${key}_${lval}
    echo "$((i++))"
    runname=post_bout_${key}_${lval}

    qsub -q serial -pe 1way 16 -N $runname -cwd -V -m e -A BOUT++_startup -l h_rt=2:00:00 -j y wrap3.py $current_dir
    
    
    allruns=$runname,$allruns
done
#echo ${allruns[@]}
allruns=${allruns/%,/} #get rid of the last comma

echo $allruns
qsub -q serial -pe 1way 16 -N "post_bout_re_eig64" -cwd -V -m e -A \
    BOUT++_startup -l h_rt=4:00:00 -j y \
    -hold_jid $allruns wrap.py #make this take an input var for the status.log 

#qsub -q serial -pe 1way 16 -N "post_bout" -cwd -V -m e -A BOUT++_startup -l h_rt=4:00:00 -j y wrap3.py
#-check the result
#idl runidl.pro
