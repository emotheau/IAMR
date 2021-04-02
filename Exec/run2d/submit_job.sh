#!/bin/bash

##SBATCH --qos=regular
##SBATCH --time=00:30:00
##SBATCH --nodes=2
##SBATCH --tasks-per-node=64
##SBATCH --constraint=haswell

res=256
chain_num_current=$1
chain_num_next=$(($1+1))
EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
OUTPUTDIR_prev='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/'$res/$chain_num_current
OUTPUTDIR_next='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/'$res/$chain_num_next

mkdir -p $OUTPUTDIR_next

INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
#RESTART='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/CHECKPOINT_FILES/NX_256/'
echo $OUTPUTDIR_prev
echo $OUTPUTDIR_next
rm -r $OUTPUTDIR_prev/*.temp
#TODO the regex is not reliable especially when the timestep iter crosses 99999. until it is fixed be careful in specifying num in [0-9]{num,}
latest_checkpoint=$(find $OUTPUTDIR_prev -maxdepth 1 -print | grep -E 'truth[0-9]{6,}' | sort -r | head -1) 
echo $latest_checkpoint

RESTART=/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/restartdir/$res
#mkdir -p $RESTART
cp -TR $latest_checkpoint $RESTART/latest_checkpoint

echo | ls $RESTART/latest_checkpoint

srun -n 128 $EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.Kolmogorov amr.check_file=$OUTPUTDIR_next/truth amr.restart=$RESTART/latest_checkpoint amr.check_per=0.1 amr.n_cell= $res $res amr.v=0 ns.v=0
