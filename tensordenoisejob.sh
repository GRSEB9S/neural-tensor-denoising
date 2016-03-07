#!/bin/bash
#$ -S /bin/bash
#$ -l mem=8G,time=24:00:00
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
echo "Your email: $jss2219@cumc.columbia.edu"

# change directory
cd /ifs/scratch/zmbbi/la_lab/jss2219/

#The following is the job to be performed:
/nfs/apps/matlab/2015a/bin/matlab -singleCompThread -nojvm -nodisplay -nosplash -r '/ifs/scratch/zmbbi/la_lab/jss2219/tensorDenoiseCluster($SGE_TASK_ID)'

