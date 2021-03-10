#!/bin/bash

basedir=/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/256

restart_dir=/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/restartdir/256

filename=$(find $basedir -maxdepth 1 -print | grep -E 'truth[0-9]{6,}' | sort -r | head -1) 
echo $filename
#mkdir -p $restart_dir
#
#cp -r $filename $restart_dir'/latest_checkpoint'
#echo $filename
