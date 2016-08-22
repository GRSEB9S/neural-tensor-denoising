#!/bin/bash
#PBS -N dat2plot
#PBS -W group_list=yetizmbbi
#PBS -l walltime=24:00:00,mem=4000mb
#PBS -M jss2219@cumc.columbia.edu
#PBS -m n
#PBS -V
#PBS -t 1-5
#PBS -o /vega/zmbbi/users/jss2219/TensorDenoising/yeti-output/
#PBS -e /vega/zmbbi/users/jss2219/TensorDenoising/yeti-error/

#The following is the job to be performed:
matlab -nosplash -nodisplay -nodesktop -r 'tensorDenoiseDat2Plot($PBS_ARRAYID)' >/vega/zmbbi/users/jss2219/TensorDenoising/mat-outfile/matplot$SGE_TASK_ID

