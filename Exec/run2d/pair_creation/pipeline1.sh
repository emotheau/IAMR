#!/bin/bash
TRUTHDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/timestep_0_1'
COARSEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/process/timestep_0_1/coarse'
TOOLDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Util/ConvertCheckpoint'
EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
TIMESTEPDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/process/timestep_0_1/timestepped'
REFINEDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/process/timestep_0_1/refined'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
PLOTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/process/timestep_0_1/plots'

filename=$1

#$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.Kolmogorov  amr.restart=$filename amr.plotfile_on_restart=1 amr.regrid_on_restart=1 amr.n_cell= 256 256 amr.plot_file="$PLOTDIR/plots"

$TOOLDIR/ConvertCheckpointGrids2d.gnu.haswell.ex checkin=$filename checkout="$COARSEDIR/coarse_${filename##*/}" interp_kind=coarsen user_ratio=2

#wait

$EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.KolmogorovPairs amr.restart="$COARSEDIR/coarse_${filename##*/}" amr.plotfile_on_restart=0 stop_interval=0.1 amr.n_cell=128 128 amr.check_file="$TIMESTEPDIR/tilde_${filename##*/}_" 

