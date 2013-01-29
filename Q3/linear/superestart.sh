
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


#llist=(.5)
#llist=(1 .8 .7 .5 .3 .1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4)
#llist=(1 .7 .8 .5 .3)
#llist=(.1 5e-2 1e-2)
#tstep=(1e-3 1e-3 1e-3 1e-3 1e-3 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2) 
#1e-3)
#1e-2 1e-2 1e-2 1e-2 1e-2)

llist=(1 .5 .1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4)

tstep=(1 1 1 1 1 1 1 1 1) 

key=gam64_detail
queue=normal


rm status_re_$key.log

make
i=0
allruns=()
for lval in ${llist[@]}
do
    current_dir=$data_dir/data_re_${key}_${lval}
    IC_dir=$data_dir/data_${key}_${lval}

    rm -r $current_dir
    mkdir $current_dir

    ls $IC_dir | grep BOUT.restart
    cp $IC_dir/BOUT.restart.*.nc $current_dir
    echo $IC_
    echo $current_dir

  #  echo "diff $IC_dir/BOUT.restart.0.nc $current_dir/BOUT.restart.0.nc"
  #  diff $IC_dir/BOUT.restart.0.nc $current_dir/BOUT.restart.0.nc
    #echo $current_dir
    runname=q3_simp_re_${key}_${lval}

    #rm -r $current_dir
    

  
  #ltime=$(echo "scale=5; ${lval}*.02" | bc -l)
  
    cp q3_simp.cxx $current_dir/q3_simp.cxx.ref
    #cp q3_simp.cxx   $PWD/data_${lval}/q3_simp.cxx.ref

    
    sed "s/ZMAX = 1/ZMAX = ${lval}/g" BOUT.inp > temp.inp
    #sed "s/ys_opt = 2/ys_opt = 3/g" temp.inp > temp2.inp
    #sed "s/xs_opt = 1/xs_opt = 2/g" temp.inp > temp2.inp 
    #sed "s/xs_mode = 2/xs_opt = 1/g" temp.inp > temp2.inp
    #sed "s/TIMESTEP = 1/TIMESTEP =  1e-2/g" temp.inp > $current_dir/BOUT.inp
    sed "s/restart = false/restart = true/g" temp.inp > temp2.inp
    sed "s/NOUT = 100/NOUT = 200/g" temp2.inp > temp3.inp

    #sed "s/ZMAX = .0001/ZMAX = ${lval}/g" BOUT.inp > temp.inp
    #sed "s/TIMESTEP = 1/TIMESTEP =  ${tstep[$i]}/g" temp2.inp > $current_dir/BOUT.inp
    sed "s/TIMESTEP = 1/TIMESTEP =  ${tstep[$i]} /g" temp3.inp > $current_dir/BOUT.inp
    echo "$((i++))"
    
    ln -s $current_dir $PWD/data_re_${key}_${lval}
    echo $current_dir >> status_re_${key}.log

    if [ "${lval}" = "${llist[${#llist[@]}-1]}" ]; then
	echo 'done' #>> status.log #for parallel submission this makes no sense
    fi
    
    #ibrun -n $NP -o 0  ./q3_simp -d $current_dir
    qsub -q $queue -pe 16way 16 -N $runname -cwd -V -m e -A BOUT++_startup -l h_rt=4:00:00 -j y ./subrun.sh $NP ./q3_simp $current_dir -re
    
    
    allruns=$runname,$allruns
done
#echo ${allruns[@]}
allruns=${allruns/%,/} #get rid of the last comma
diff $IC_dir/BOUT.restart.0.nc $current_dir/BOUT.restart.0.nc
echo $allruns
qsub -q serial -pe 1way 16 -N "post_bout_gam64" -cwd -V -m e -A \
    BOUT++_startup -l h_rt=4:00:00 -j y \
    -hold_jid $allruns wrap.py status_${key}.log output_$key.pdf

 #make this take an input var for the status.log 

#qsub -q serial -pe 1way 16 -N "post_bout" -cwd -V -m e -A BOUT++_startup -l h_rt=4:00:00 -j y wrap.py
#-check the result
#idl runidl.pro
