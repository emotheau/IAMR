#!/bin/bash

res=256
iname="hit_ic_4_"$res".dat"

basedir="/global/cscratch1/sd/jpathak/iamr_hit/"
#basedir="/project/projectdirs/dasrepo/jpathak/IAMR-HIT/Exec/HIT/"
mkdir -p $basedir

check_path=$basedir$res"/checkpoints/"
check_file=$check_path"/chk_"
plot_path=$basedir$res"/plotfiles/"
plot_file=$plot_path"/plt_"

verbose=1


srun -n 512 amr3d.gnu.haswell.MPI.ex inputs.3d.HIT prob.inres=$res prob.iname=$iname amr.n_cell=$res $res $res amr.plot_file=$plot_file amr.check_file=$check_file amr.v=$verbose ns.v=$verbose

