#!/bin/bash
#$ -S /bin/bash
#$ -l mem=4G,time=12:00:00
#$ -cwd
#$ -N tensorDenoising
#$ -j y
#$ -m ea # send email at end of job and if aborted
#$ -t 1-240

## Some Debugging info
echo "Your UNI: $USER"
echo “Starting on : $(date)”
echo “Running on node : $(hostname)”
echo “Current directory : $(pwd)”
echo “Current job ID : $JOB_ID”
echo “Current job name : $JOB_NAME”
# echo "Your email: $USER@c2b2.columbia.edu"

#The following is the job to be performed:
/nfs/apps/matlab/2015a/bin/matlab -singleCompThread -nojvm -nodisplay -nosplash -r 'tensorDenoiseCluster($SGE_TASK_ID)' >matoutfile$SGE_TASK_ID

