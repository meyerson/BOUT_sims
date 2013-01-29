
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
make

#-run the case       
echo Running with NP = $NP       

current_dir=$PWD
echo $current_dir


data_dir=`echo $PWD|sed "s/share\/home/scratch/g"`

NOUTS=(1000 1000 1000 200)
tstep=(1e0 1e-2 1e2)
#llist=(5e-2 1e-1 5e-3)
llist=(1e0)
#second set of jobs, submitted to development q
tstep2=(5e-1)
llist2=(8e-2)
queue=development

old_key=2dfast
key=2dfast_re
rm status_$key.log

make
i=0
allruns=()
for lval in ${llist[@]}
do
    current_dir=$data_dir/data_${key}_${lval}
    echo $current_dir
    runname=drift_${key}_${lval}


    
    #rm -r $current_dir
    mkdir -p $current_dir

    
    cp hlmk.cxx $current_dir/physics_code.cxx.ref
    cp $data_dir/data_${old_key}_${lval}/*restart* $current_dir/

#cp 2fluid.cxx   $PWD/data_${lval}/2fluid.cxx.ref

    sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT.inp > temp.inp
    sed "s/NOUT = 100/NOUT = ${NOUTS[$i]}/g" temp.inp > temp2.inp
    sed "s/restart = false/restart = true/g" temp2.inp > temp3.inp
    sed "s/nonlinear = false/nonlinear = true/g" temp3.inp > temp4.inp
    #sed "s/par_damp = true/par_damp = false/g" temp4.inp > temp5.inp

    sed "s/TIMESTEP = 1/TIMESTEP =  ${tstep[$i]}/g" temp4.inp > $current_dir/BOUT.inp
    
    echo "$((i++))"
    
    ln -s $current_dir $PWD/data_${key}_${lval}
    echo $current_dir >> status_${key}.log

    if [ "${lval}" = "${llist[${#llist[@]}-1]}" ]; then
	echo 'done' #>> status.log #for parallel submission this makes no sense
    fi
    
    #ibrun -n $NP -o 0  ./2fluid -d $current_dir
    qsub -q $queue -pe 16way 16 -N $runname -cwd -V -m e -A BOUT++_startup -l h_rt=2:00:00 -j y ./subrun.sh $NP ./hlmk $current_dir re
    
    
    allruns=$runname,$allruns
done
#echo ${allruns[@]}
allruns=${allruns/%,/} #get rid of the last comma


 qsub -q serial -pe 1way 16 -N "post_bout"_$key -cwd -V -m e -A \
     BOUT++_startup -l h_rt=4:00:00 -j y \
     -hold_jid $allruns wrap.py status_${key}.log output_$key.pdf
 #make this take an input var for the status.log 

#qsub -q serial -pe 1way 16 -N "post_bout"_$key  -cwd -V -m e -A BOUT++_startup -l h_rt=4:00:00 -j y wrap.py status_${key}.log output_$key.pdf
