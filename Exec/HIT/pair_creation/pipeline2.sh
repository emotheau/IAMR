#!/bin/bash


idx=$1

res=256
ur=4

TOOLDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Util/ConvertCheckpoint'
EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
TIMESTEPDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'$res'/ur'$ur'/'$idx'/timestepped'
REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'$res'/ur'$ur'/'$idx'/refined'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
REFINEDPLOTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'$res'/ur'$ur'/'$idx'/refined_plots'


filename=$2

$TOOLDIR/ConvertCheckpointGrids2d.gnu.haswell.ex checkin=$filename checkout="$REFINEDIR/refined_${filename##*/}" interp_kind=refine user_ratio=$ur

$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPlot  amr.restart="$REFINEDIR/refined_${filename##*/}" amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= $res $res amr.plot_file="$REFINEDPLOTDIR/plots_${filename##*/}" amr.max_grid_size=$res amr.blocking_factor=$res

