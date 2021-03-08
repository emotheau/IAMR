#!/bin/bash

EXECDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
OUTPUTDIR='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/timestep_0_1'
INPUTDIR='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'
RESTART='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d'

mkdir -p $OUTPUTDIR

srun -n 256 $EXECDIR/amr2d.gnu.haswell.MPI.ex $INPUTDIR/inputs.2d.Kolmogorov amr.check_file=$OUTPUTDIR/truth amr.restart=$RESTART/chk20000

