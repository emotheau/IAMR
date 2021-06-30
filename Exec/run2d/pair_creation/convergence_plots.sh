#!/bin/bash

idx=$1

res=1024
ur=16
l_res=$(($res/$ur))

EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
PLOTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/convergence_test/'$res'/'

filename='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/CHECKPOINT_FILES/NX_1024/chk144566/'

$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPlot  amr.restart=$filename amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= $res $res amr.plot_file="$PLOTDIR/plots" amr.max_grid_size= $res amr.blocking_factor= $res

#$TOOLDIR/ConvertCheckpointGrids2d.gnu.haswell.ex checkin=$filename checkout="$COARSEDIR/coarse_${filename##*/}" interp_kind=coarsen user_ratio=$ur

#wait

#$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPairs amr.restart="$COARSEDIR/coarse_${filename##*/}" amr.plotfile_on_restart=-1 stop_interval=0.1 amr.n_cell=$l_res $l_res amr.check_file="$TIMESTEPDIR/tilde_${filename##*/}_" amr.plot_int=10000000 amr.plot_file="plotfiles/" 

