#!/usr/bin/sh

#module purge 
#module load cuda/11.2.1-gcc-10.2.0

#source /projects/b1021/Jianshe/codes/anaconda3/etc/profile.d/conda.sh
#conda activate hoomd-env3_13

# alternatively, use sys,path.append in the script to point to the MPI or non-MPI version

#export PYTHONPATH=$PYTHONPATH:/projects/b1021/Jianshe/codes/hoomd-v3_13/install

fname=create_bs_topology.py
outfname=branch_init

python $fname $outfname
python modify_gsd.py $outfname $outfname
