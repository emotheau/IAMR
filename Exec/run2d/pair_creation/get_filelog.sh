#!/bin/bash


CHKDIR=$1   #'/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/timestep_0_1'

for filename in $CHKDIR/*; do
  echo $filename >> files.txt
done
