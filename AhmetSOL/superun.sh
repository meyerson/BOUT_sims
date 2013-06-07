
#!/bin/bash

NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    NP=132
    #NP=48
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

current_dir=$PWD
echo $current_dir
queue=normal


data_dir=`echo $PWD|sed "s/home1/scratch/g"`

NOUTS=(350 510)
tstep=(1e1)
#llist=(1e-2)
llist=(2e-3)

key=AhmetInBout_HD
rm status_$key.log

make
i=0
allruns=()
for lval in ${llist[@]}
do
    current_dir=$data_dir/data_${key}_${lval}_HD
    echo $current_dir
    runname=${key}_${lval}_${NP}

    #rm -r $current_dir
    mkdir -p $current_dir

    
    cp SOLblob.cxx $current_dir/physics_code.cxx.ref
   
    sed "s/NOUT = 100/NOUT = ${NOUTS[$i]}/g" BOUT_${key}.inp > temp.inp
    #sed "s/NXPE = 264/NXPE = ${NP}/g" temp.inp > temp2.inp
    sed "s/mu = 1e-1/mu = ${lval}/g" temp.inp > temp2.inp
    sed "s/nu = 1e-1/nu = ${lval}/g" temp2.inp > temp.inp
    sed "s/NXPE = 1/NXPE = ${NP}/g" temp.inp > temp2.inp

    sed "s/TIMESTEP = 1/TIMESTEP =  ${tstep[$i]}/g" temp2.inp > $current_dir/BOUT.inp
    
    echo "$((i++))"
    
    ln -s $current_dir $PWD/data_${key}_${lval}
    echo $current_dir >> status_${key}.log

    if [ "${lval}" = "${llist[${#llist[@]}-1]}" ]; then
	echo 'done' #>> status.log #for parallel submission this makes no sense
    fi
    
    #ibrun -n $NP -o 0  ./2fluid -d $current_dir
    qsub -q $queue -pe 12way $NP -N $runname -cwd -V -m e -A BOUT++_startup -l h_rt=5:00:00 -j y ./subrun.sh $NP ./SOLblob $current_dir
    
    
    allruns=$runname,$allruns
done
#echo ${allruns[@]}
allruns=${allruns/%,/} #get rid of the last comma


 #qsub -q serial -pe 1way 16 -N "post_bout"_$key -cwd -V -m e -A \
 #    BOUT++_startup -l h_rt=4:00:00 -j y \
 #    -hold_jid $allruns wrap.py status_${key}.log output_$key.pdf
 #make this take an input var for the status.log 

#qsub -q serial -pe 1way 16 -N "post_bout"_$key  -cwd -V -m e -A BOUT++_startup -l h_rt=4:00:00 -j y wrap.py status_${key}.log output_$key.pdf
