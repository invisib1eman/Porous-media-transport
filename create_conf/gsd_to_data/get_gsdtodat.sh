#!/bin/bash

module purge all
source /projects/b1021/Jianshe/codes/anaconda3/etc/profile.d/conda.sh
conda activate hoomd-env3_13

export pre_gsd="branch_init"
export mdf_gsd="temp"
export data_file="branch_rba"

#python modify_gsd.py $pre_gsd.gsd $mdf_gsd.gsd

python gsdtodat.py $pre_gsd.gsd $data_file.data

#conda deactivate
#conda activate ovito-env

#python outData_ovito.py $mdf_gsd.gsd $data_file.data

#rm $mdf_gsd.gsd
