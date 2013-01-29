
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
llist=(1 .5 .1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4)
#.5 .1 .2 .3)
#llist=(.05 
#llist=(.5 .3 .2 .004) 
#llist = (8e-4 5e-4)
#llist=(.3 5e-4)
tstep=(1 1 1 1e1 1e1 1e1 1e1 1e1 1e1) 
#1e-3)
#1e-2 1e-2 1e-2 1e-2 1e-2)
key=gam64_detail
queue=normal

#key=fast_g64
rm status_$key.log

make
i=0
allruns=()
for lval in ${llist[@]}
do
    current_dir=$data_dir/data_${key}_${lval}
    echo $current_dir
    runname=q3_simp_${key}_${lval}

    #rm -r $current_dir
    mkdir $current_dir

    echo $ltime
    
    cp q3_simp.cxx $current_dir/q3_simp.cxx.ref
    #cp q3_simp.cxx   $PWD/data_${lval}/q3_simp.cxx.ref

    
    sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT.inp > temp.inp
    #sed "s/ys_opt = 2/ys_opt = 3/g" temp.inp > temp2.inp
    #sed "s/xs_opt = 1/xs_opt = 2/g" temp.inp > temp2.inp 
    #sed "s/xs_mode = 2/xs_opt = 1/g" temp.inp > temp2.inp
    #sed "s/TIMESTEP = 1/TIMESTEP =  ${ltime}/g" temp.inp > $current_dir/BOUT.inp

    #sed "s/ZMAX = .0001/ZMAX = ${lval}/g" BOUT.inp > temp.inp
    #sed "s/TIMESTEP = 1/TIMESTEP =  ${tstep[$i]}/g" temp2.inp > $current_dir/BOUT.inp
    sed "s/TIMESTEP = 1/TIMESTEP =  ${tstep[$i]}/g" temp.inp > $current_dir/BOUT.inp
    echo "$((i++))"
    
    ln -s $current_dir $PWD/data_${key}_${lval}
    echo $current_dir >> status_${key}.log

    if [ "${lval}" = "${llist[${#llist[@]}-1]}" ]; then
	echo 'done' #>> status.log #for parallel submission this makes no sense
    fi
    
    #ibrun -n $NP -o 0  ./q3_simp -d $current_dir
    qsub -q $queue -pe 16way 16 -N $runname -cwd -V -m e -A BOUT++_startup -l h_rt=2:00:00 -j y ./subrun.sh $NP ./q3_simp $current_dir
    
    
    allruns=$runname,$allruns
done
#echo ${allruns[@]}
allruns=${allruns/%,/} #get rid of the last comma

echo $allruns
 qsub -q serial -pe 1way 16 -N "post_bout"_$key -cwd -V -m e -A \
     BOUT++_startup -l h_rt=4:00:00 -j y \
     -hold_jid $allruns wrap.py status_${key}.log output_$key.pdf
 #make this take an input var for the status.log 

#qsub -q serial -pe 1way 16 -N "post_bout"_$key  -cwd -V -m e -A BOUT++_startup -l h_rt=4:00:00 -j y wrap.py status_${key}.log output_$key.pdf
#-check the result
#idl runidl.pro
