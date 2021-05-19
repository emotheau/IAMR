#!/bin/bash


module load parallel


#idx=3
res=256
ur=4
for idx in 8 9 10 11 12 13 14
do
  rm logfile_pipe1.txt
  rm logfile_pipe2.txt
  rm files.txt
  
  TRUTHDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/'$res'/'$idx
  REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'$res'/ur'$ur'/'$idx'/refined'
  TIMESTEPDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'$res'/ur'$ur'/'$idx'/timestepped'
  
  ./get_filelog.sh $TRUTHDIR
  
  #wait
  #run gnuparallel on first part of pipeline
  #
  srun parallel --jobs 50 --resume-failed --joblog logfile_pipe1.txt ./pipeline1.sh $idx {} < files.txt
  
#  wait
  
  rm files.txt
  
#  wait
  /get_filelog.sh $TIMESTEPDIR
  
  #wait
  #run gnuparallel on second part of pipeline
  
  srun parallel --jobs 50 --resume-failed --joblog logfile_pipe2.txt ./pipeline2.sh $idx {} < files.txt
  
#  rm logfile_pipe3.txt
#  #option to plot the unrefined timestepped checkpoints
#  srun parallel --jobs 50 --resume-failed --joblog logfile_pipe3.txt ./pipeline3.sh $idx {} < files.txt
done
