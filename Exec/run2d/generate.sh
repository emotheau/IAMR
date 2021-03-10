#!/bin/bash

EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
OUTPUTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/2048'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
#RESTART='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/CHECKPOINT_FILES/NX_256/'

mkdir -p $OUTPUTDIR

latest_checkpoint=$(find $OUTPUTDIR -maxdepth 1 -print | grep -E 'truth[0-9]{6,}' | sort -r | head -1) 
RESTART=/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/restartdir/2048
mkdir -p $RESTART
cp -r $latest_checkpoint $RESTART'/latest_checkpoint'

echo $latest_checkpoint
echo $RESTART/latest_checkpoint
echo $(ls $RESTART/latest_checkpoint) 

srun -n 256 $EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.Kolmogorov amr.check_file=$OUTPUTDIR/truth amr.restart=$RESTART/latest_checkpoint amr.check_per=0.1 amr.n_cell= 2048 2048
