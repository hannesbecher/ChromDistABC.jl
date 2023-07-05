#!/bin/bash

# Grid Engine options
#$ -cwd
#$ -t 1-1000
#$ -l h_rt=05:00:00
#$ -l h_vmem=40G
#$ -N abcMouseCount10x1000               


id=`printf %04d $SGE_TASK_ID`
echo "==========================================================="
echo Running task $id on $HOSTNAME
echo "==========================================================="

SECONDS=0
bash juliaAbcCountsArg.sh human.csv humanCounts$id 10 6
duration=$SECONDS
echo "$(($duration/60))m$(($duration%60))s elapsed"
echo "==========================================================="

