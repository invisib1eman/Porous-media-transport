#!/bin/bash

#SBATCH --job-name="GN-R0.8"
#SBATCH -A b1021 	               
#SBATCH -p buyin  	 
#SBATCH -t 48:00:00               
#SBATCH -N 1 		                  
#SBATCH --ntasks-per-node=24     
#SBATCH -n 24          
#SBATCH -e errlog.%j 
#SBATCH -o outlog.%j  

source /projects/b1021/Jianshe/codes/lammps/lammps-stable/install/bin/lammps-29Sep2021.sh

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=1

dirs=("CreateNET" "EquilGel" "EquilNet")
lammps_script=("in.network_create" "in.gel_equil" "in.network_equil")

execuable="lmp_exchange_cpu"
seed=`shuf -i 100000-999999 -n 1`
init_data="branch_rba.data"

for (( i=0; i<2; i++ ))
do
  if [ $i -eq 0 ]
  then 
    mkdir ${dirs[i]}
    mv ${lammps_script[i]} ${init_data} ${dirs[i]}/
    cd ${dirs[i]}/
    mpirun ${execuable} -v seed ${seed} -in ${lammps_script[i]} 
    wait
    cd ..
  elif [ $i -eq 1 ]
  then
    mkdir ${dirs[i]}/
    mv ${lammps_script[i]} ${dirs[i]}/
    cp ${dirs[0]}/network.data ${dirs[i]}/
    cd ${dirs[i]}/
    mpirun ${execuable} -v seed ${seed} -in ${lammps_script[i]} 
    wait
    cd ../
  fi
done






