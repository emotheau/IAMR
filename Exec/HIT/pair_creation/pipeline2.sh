#!/bin/bash

res=128
ur=2
l_res=$(($res/$ur))

BASEDIR='/project/projectdirs/dasrepo/jpathak/IAMR-HIT/Exec/HIT/'
DATADIR='/global/cscratch1/sd/jpathak/iamr_hit/'$res

TRUTHDIR=$DATADIR'/checkpoints'
COARSEDIR=$DATADIR'/coarse'
TIMESTEPDIR=$DATADIR'/timestepped'
REFINEDIR=$DATADIR'/refined'
REFINEDPLOTDIR=$DATADIR'/refined_plot'
PLOTDIR=$DATADIR'/plots'

INPUTDIR=$BASEDIR
TOOLDIR=$BASEDIR
EXECDIR=$BASEDIR


filename='/global/cscratch1/sd/jpathak/iamr_hit/128/timestepped/tilde_chk_00475_00480/'

#$TOOLDIR/ConvertCheckpointGrids3d.gnu.haswell.ex checkin=$filename checkout="$REFINEDIR/refined_${filename##*/}" interp_kind=refine user_ratio=$ur

$EXECDIR/amr3d.gnu.haswell.MPI.ex $INPUTDIR/inputs.3d.HIT  amr.restart="$REFINEDIR/refined_${filename##*/}" amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= $res $res $res amr.plot_file="$REFINEDPLOTDIR/plots_${filename##*/}" amr.max_grid_size=$res amr.blocking_factor=$res

