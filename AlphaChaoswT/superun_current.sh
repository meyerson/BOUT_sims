#!/bin/bash
NO_ARGS=0 
OPTERROR=65

if [ $# -eq "$NO_ARGS" ]  # Script invoked with no command-line args?
then
    #NP=108
    #NP=48
    NP=72
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
#queue=normal
queue=development
#project=A-ph5
project=BOUT++_startup

data_dir=`echo $PWD|sed "s/home1/scratch/g"`

NOUTS=(3000 3000 3000 3000 3000)
NOUT=100
# 350 350 350)
tsteps=(3e0 3e0 3e0 1e0 1e0)
tstep=3e-1

alphalist=(true true true)
solverlist=(cvode cvode cvode cvode cvode)
solver=cvode
#.04 ~ 30
mulist=(1e-2 1e-2 1e-2 1e-2 1e-2)
nulist=(1e-2 1e-2 1e-2 1e-2 1e-2)
mu=1e-2
nu=1e-2
chaosalpha=(chaos smooth)

alphas=(3.0e-3 3.0e-3)
alpha=3.0e-3

betas=(1.0e-2 1.0e-2 1.0e-2 1.0e-2 1.0e-2 1.0e-2)
beta=1.0e-2


eps=.1
#dx=.267
dx=.1
#nx=2164
nx=1084
##eps=1.0e-1
#dx=.2
ZMAX=10.0

dx_list=(.04 .04 .2 .2 .2 .2 .2)

nlog=true
#(true true true false false false)
sol_edge=0.0
inc_jpar=false
inc_jpar_list=(true true true false false false)
#rho = .4 cm
key=SOLwT_dev
rm status_$key.log

make
allruns=()


i=0

for a_map in ${chaosalpha[@]}
do

    current_dir=$data_dir/data_${key}_nlog_${a_map}_a${alpha}_eps${eps}_dev_dx${dx}

    echo $current_dir
    #runname=SOL_${chaosalpha[$i]}_${alpha}
    runname=${key}_nlog_${alpha}_${NP}
    #rm -r $current_dir
    mkdir -p $current_dir

    cp SOLblob.cxx $current_dir/physics_code.cxx.ref
   
    sed "s/NOUT = 100/NOUT = ${NOUT}/g" BOUT_${key}.inp > temp.inp
  
    sed "s/mu = 1e-1/mu = ${mu}/g" temp.inp > temp2.inp
    sed "s/nu = 1e-1/nu = ${nu}/g" temp2.inp > temp.inp
    sed "s/NXPE = 1/NXPE = ${NP}/g" temp.inp > temp2.inp
    
    sed "s/chaosalpha = smooth/chaosalpha = ${a_map}/g" temp2.inp > temp.inp

    sed "s/alpha_c = 1.0e-3/alpha_c = ${alpha}/g" temp.inp > temp2.inp
    sed "s/TIMESTEP = 5e0/TIMESTEP =  ${tstep}/g" temp2.inp > temp.inp
    sed "s/dx = .1/dx = ${dx}/g" temp.inp > temp2.inp
    sed "s/beta = 2.0e-2/beta = ${beta}/g" temp2.inp > temp.inp
    sed "s/eps = .1/eps = ${eps}/g" temp.inp > temp2.inp
    sed "s/log_n = true/log_n = ${nlog}/g" temp2.inp > temp.inp
    
    #sed "s/inc_jpar = false/inc_jpar = ${inc_jpar}/g" temp2.inp > temp.inp
    sed "s/type = cvode/type = ${solver}/g" temp.inp > temp2.inp
    sed "s/ZMAX = 1/ZMAX = ${ZMAX}/g" temp2.inp > temp.inp
    sed "s/nx = 724/nx = ${nx}/g" temp.inp > $current_dir/BOUT.inp
    
    echo "$((i++))"
    
    #ln -s $current_dir $PWD/data_${key}_${lval}
    echo $current_dir >> status_${key}.log

    qsub -q $queue -pe 12way ${NP} -N $runname -cwd -V -m e -A $project -l h_rt=1:00:00 -j y ./subrun.sh $NP ./SOLblobwT $current_dir
    
    
    allruns=$runname,$allruns
done

allruns=${allruns/%,/} #get rid of the last comma

