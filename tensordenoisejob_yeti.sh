#!/bin/bash
#PBS -N tensordenoising
#PBS -W group_list=yetizmbbi
#PBS -l walltime=24:00:00,mem=4000mb
#PBS -M jss2219@cumc.columbia.edu
#PBS -m n
#PBS -V
#PBS -t 1-500
#PBS -o /vega/zmbbi/users/jss2219/TensorDenoising/yeti-output/
#PBS -e /vega/zmbbi/users/jss2219/TensorDenoising/yeti-error/


#Apparently matlab only works in an interactive session.
qsub -q interactive -I -W group_list=yetizmbbi -l walltime=04:00:00,mem=4000mb

#Run Matlab
matlab -nosplash -nodisplay -nodesktop -r "tensorDenoiseCluster($PBS_ARRAYID)" > /vega/zmbbi/users/jss2219/TensorDenoising/mat-outfile/matoutfile$PBS_ARRAYID
