#!/bin/bash
TRUTHDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/2048'
COARSEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/2048/coarse'
TOOLDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Util/ConvertCheckpoint'
EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
TIMESTEPDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/2048/timestepped'
REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/2048/refined'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
PLOTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/test/2048/'

filename=$TRUTHDIR/truth281794
#filename='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/CHECKPOINT_FILES/NX_2048/chk281681'
$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.Kolmogorov  amr.restart=$filename amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= 2048 2048 amr.plot_file="$PLOTDIR/plots" amr.max_grid_size= 2048 amr.blocking_factor= 2048 max_step=0


#$TOOLDIR/ConvertCheckpointGrids2d.gnu.haswell.ex checkin=$filename checkout="$COARSEDIR/coarse_${filename##*/}" interp_kind=coarsen user_ratio=4
#
##wait
#
#$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPairs amr.restart="$COARSEDIR/coarse_${filename##*/}" amr.plotfile_on_restart=0 stop_interval=0.1 amr.n_cell=512 512 amr.check_file="$TIMESTEPDIR/tilde_${filename##*/}_" 
#
