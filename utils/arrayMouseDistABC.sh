#!/bin/bash

# Grid Engine options
#$ -cwd
#$ -t 1-1000
#$ -l h_rt=01:00:00
#$ -l h_vmem=10G
#$ -N abcMouseDist10x1000 


id=`printf %04d $SGE_TASK_ID`
echo "==========================================================="
echo Running task $id on $HOSTNAME
echo "==========================================================="

SECONDS=0
bash juliaAbcDistancesArg.sh mouse.csv NEWmouseDist$id 100 8
duration=$SECONDS
echo "$(($duration/60))m$(($duration%60))s elapsed"
echo "==========================================================="

