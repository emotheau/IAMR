#!/bin/bash

TOOLDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Util/ConvertCheckpoint'
EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
TIMESTEPDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/process/timestep_0_1/timestepped'
REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/process/timestep_0_1/refined'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
REFINEDPLOTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/process/timestep_0_1/refined_plots'


filename=$1

$TOOLDIR/ConvertCheckpointGrids2d.gnu.haswell.ex checkin=$filename checkout="$REFINEDIR/refined_${filename##*/}" interp_kind=refine user_ratio=2

$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPair  amr.restart="$REFINEDIR/refined_${filename##*/}" amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= 256 256 amr.plot_file="$REFINEDPLOTDIR/plots_${filename##*/}"

