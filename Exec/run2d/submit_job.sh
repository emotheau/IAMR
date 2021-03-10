#!/bin/bash

#SBATCH --qos=regular
#SBATCH --time=00:03:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell

res=512

EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
OUTPUTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/'$res
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
#RESTART='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/CHECKPOINT_FILES/NX_256/'

rm -r $OUTPUTDIR/*.temp
latest_checkpoint=$(find $OUTPUTDIR -maxdepth 1 -print | grep -E 'truth[0-9]{5,}' | sort -r | head -1) 
echo $latest_checkpoint
RESTART=/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/restartdir/$res
mkdir -p $RESTART
cp -TR $latest_checkpoint $RESTART/latest_checkpoint

echo $latest_checkpoint
echo $RESTART/latest_checkpoint
echo $(ls $RESTART/latest_checkpoint) 

srun -n 32 $EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.Kolmogorov amr.check_file=$OUTPUTDIR/truth amr.restart=$RESTART/latest_checkpoint amr.check_per=0.1 amr.n_cell= $res $res amr.v=0 ns.v=0
