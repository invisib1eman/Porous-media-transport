#!/bin/bash

#SBATCH --job-name="equilGEL"
#SBATCH -A b1030 	               
#SBATCH -p buyin  	 
#SBATCH -t 48:00:00               
#SBATCH -N 1 		                  
#SBATCH --ntasks-per-node=28     
#SBATCH -n 28          
#SBATCH -e errlog.%j 
#SBATCH -o outlog.%j  

source /projects/b1021/Jianshe/codes/lammps/lammps-stable/install/bin/lammps-29Sep2021.sh

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=1

resume=1
gpu_flag=0
data_flag=1 # 0 for data, 1 for restart
npt_flag=0 # run npt
equil_flag=0 # equilibate the system

nsteps_output=400000000
nsteps_equil=10000000

rba=0.67
energy_barrier=0.0
temp=1.0
den=1.1
press=12

jobname="gel_equil"
simname="branch_rba"${rba}"_p"${press}"_T"${temp}
input_fname="2.data"
seed=`shuf -i 100000-999999 -n 1`

if [ $gpu_flag -eq 1 ] 
then
    execuable="lmp_exchange"
else
    execuable="lmp_exchange_cpu"
fi


if [ $resume -eq 0 ]
then
    if [ $gpu_flag -eq 1 ]
    then 
        mpirun ${execuable} -pk gpu 2 -sf gpu \
        -v data_flag ${data_flag} \
        -v npt_flag ${npt_flag} \
        -v equil_flag ${equil_flag} \
        -v temp ${temp} -v press ${press} \
        -v energy_barrier ${energy_barrier} \
        -v seed ${seed} -v simname ${simname} \
        -v input_fname ${input_fname} \
        -v nsteps_equil ${nsteps_equil} -v nsteps_output ${nsteps_output} \
        -in in.${jobname}
    else
        mpirun ${execuable} \
        -v data_flag ${data_flag} \
        -v npt_flag ${npt_flag} \
        -v equil_flag ${equil_flag} \
        -v temp ${temp} -v press ${press} \
        -v energy_barrier ${energy_barrier} \
        -v seed ${seed} -v simname ${simname} -v input_fname ${input_fname} \
        -v nsteps_equil ${nsteps_equil} -v nsteps_output ${nsteps_output} \
        -in in.${jobname}
    fi
else
    if [ -e log.lammps ]
    then 
        last=`tail -n1 log.lammps | awk '{print $1}'`
        cp log.lammps log.$last
        if [ -e tmpB.sw.msd ]
        then 
            cp tmpB.sw.msd tmpB.sw.msd.$last
        fi
        resteps=$last
        nsteps_output=$(echo | awk "{print ${nsteps_output}-${resteps}}")
        t1=`stat -c %Y 1.restart`
        t2=`stat -c %Y 2.restart`
        if [ $t1 -gt $t2 ]
        then
            restartfile="1.restart"
        else
            restartfile="2.restart"
        fi
    fi

    if [ -e $restartfile ]
    then
        if [ $gpu_flag -eq 1 ]
        then 
            mpirun ${execuable} -pk gpu 2 -sf gpu \
            -v data_flag 1 \
            -v npt_flag ${npt_flag} \
            -v equil_flag ${equil_flag} \
            -v resteps ${resteps} \
            -v temp ${temp} -v press ${press} \
            -v energy_barrier ${energy_barrier} \
            -v seed ${seed} \
            -v simname ${simname} -v input_fname ${restartfile} \
            -v nsteps_equil ${nsteps_equil} -v nsteps_output ${nsteps_output} \
            -in in.${jobname}
        else
            mpirun ${execuable} \
            -v data_flag 1 \
            -v npt_flag ${npt_flag} \
            -v equil_flag ${equil_flag} \
            -v resteps ${resteps} \
            -v temp ${temp} -v press ${press} \
            -v energy_barrier ${energy_barrier} \
            -v seed ${seed} \
            -v simname ${simname} -v input_fname ${restartfile} \
            -v nsteps_equil ${nsteps_equil} -v nsteps_output ${nsteps_output} \
            -in in.${jobname}
        fi
    fi
fi
