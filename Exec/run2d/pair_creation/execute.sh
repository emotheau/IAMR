#!/bin/bash


module load parallel

rm logfile_pipe1.txt
rm logfile_pipe2.txt
rm files.txt

TRUTHDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/2048'
REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/2048/refined'
TIMESTEPDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/2048/timestepped'

./get_filelog.sh $TRUTHDIR

#wait
#run gnuparallel on first part of pipeline
#
srun parallel --jobs 32 --resume-failed --joblog logfile_pipe1.txt ./pipeline1.sh {} < files.txt

wait

rm files.txt

wait
./get_filelog.sh $TIMESTEPDIR

#wait
#run gnuparallel on second part of pipeline

srun parallel --jobs 32 --resume-failed --joblog logfile_pipe2.txt ./pipeline2.sh {} < files.txt
