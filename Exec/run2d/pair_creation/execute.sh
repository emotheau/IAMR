#!/bin/bash


module load parallel

rm logfile_pipe1.txt
rm files.txt

TRUTHDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/timestep_0_1'
REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/process/timestep_0_1/refined'

#get filelog
./get_filelog.sh $TRUTHDIR

wait
#run gnuparallel on first part of pipeline

srun parallel --jobs 50 --resume-failed --joblog logfile_pipe1.txt ./pipeline1.sh {} < files.txt

wait

rm files.txt

wait
./get_filelog.sh $REFINEDIR

wait
#run gnuparallel on second part of pipeline

srun parallel --jobs 50 --resume-failed --joblog logfile_pipe2.txt ./pipeline2.sh {} < files.txt
