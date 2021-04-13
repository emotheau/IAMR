#!/bin/bash

idx=$1

TRUTHDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/256/'$idx
COARSEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/256/ur8/'$idx'/coarse'
TOOLDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Util/ConvertCheckpoint'
EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
TIMESTEPDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/256/ur8/'$idx'/timestepped'
REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/256/ur8/'$idx'/refined'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
PLOTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/256/ur8/'$idx'/plots'

filename=$2

$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPlot  amr.restart=$filename amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= 256 256 amr.plot_file="$PLOTDIR/plots" amr.max_grid_size= 256 amr.blocking_factor= 256


$TOOLDIR/ConvertCheckpointGrids2d.gnu.haswell.ex checkin=$filename checkout="$COARSEDIR/coarse_${filename##*/}" interp_kind=coarsen user_ratio=8

#wait

$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPairs amr.restart="$COARSEDIR/coarse_${filename##*/}" amr.plotfile_on_restart=0 stop_interval=0.1 amr.n_cell=32 32 amr.check_file="$TIMESTEPDIR/tilde_${filename##*/}_" amr.plot_int=10000000 amr.plot_file="plotfiles/" 

