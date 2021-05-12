#!/bin/bash


idx=$1
res=256
ur=4

l_res=$(($res/$ur))
echo $l_res

TOOLDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Util/ConvertCheckpoint'
EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
TIMESTEPDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'$res'/ur'$ur'/'$idx'/timestepped'
REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'$res'/ur'$ur'/'$idx'/refined'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
TILDEPLOTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'$res'/ur'$ur'/'$idx'/tilde_plots'

echo $TIMESTEPDIR
echo $REFINEDIR
echo $TILDEPLOTDIR

filename=$2


$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPlot  amr.restart=$filename amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= $l_res $l_res amr.plot_file="$TILDEPLOTDIR/tilde_plots_${filename##*/}" amr.max_grid_size=$l_res amr.blocking_factor=$l_res

