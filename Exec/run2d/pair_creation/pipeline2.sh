#!/bin/bash

TOOLDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Util/ConvertCheckpoint'
EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
TIMESTEPDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/2048/timestepped'
REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/2048/refined'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
REFINEDPLOTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/2048/refined_plots'


filename=$1

#$TOOLDIR/ConvertCheckpointGrids2d.gnu.haswell.ex checkin=$filename checkout="$REFINEDIR/refined_${filename##*/}" interp_kind=refine user_ratio=4

$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPlot  amr.restart="$REFINEDIR/refined_${filename##*/}" amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= 2048 2048 amr.plot_file="$REFINEDPLOTDIR/plots_${filename##*/}" amr.max_grid_size=2048 amr.blocking_factor=2048

