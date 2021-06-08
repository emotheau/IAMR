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
PLOTDIR=$DATADIR'/plots'

INPUTDIR=$BASEDIR
TOOLDIR=$BASEDIR
EXECDIR=$BASEDIR

filename='/global/cscratch1/sd/jpathak/iamr_hit/128/checkpoints/chk_00482'
echo $filename


#srun -n 64 $EXECDIR/amr3d.gnu.haswell.MPI.ex $INPUTDIR/inputs.3d.HIT amr.restart=$filename amr.plotfile_on_restart=-1 stop_interval=0.1 amr.check_per=0.1 amr.n_cell=$res $res $res amr.check_file="$TRUTHDIR/truth" amr.plot_int=10000000 amr.plot_file="plotfiles/" 
#

srun -n 64 $EXECDIR/amr3d.gnu.haswell.MPI.ex $INPUTDIR/inputs.3d.HIT max_step=0 amr.restart=$filename amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= $res $res $res amr.plot_file="$PLOTDIR/plots" amr.max_grid_size= $res amr.blocking_factor= $res
