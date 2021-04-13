#!/bin/bash

. env_script.sh

config_name=kolmogorov256ur8crop128
run_num=1
ckpt_name=$config_name$run_num
ml_ckpt_path='/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d/ml_models/'$ckpt_name'.pt'
inference_dir='/project/projectdirs/dasrepo/jpathak/iamr_expts/inference_expts/'$config_name'_'$run_num

res=256
inp_scale_factor=0.4
target_scale_factor=4.0
ml_refinement_ratio=8
stop_interval=50

restart_ckpt=/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/256/ur8/3/coarse/coarse_truth492009/
mkdir -p $inference_dir

echo $ml_ckpt_path
echo $inference_dir
#run ml
./amr2d.gnu.haswell.MPI.ex inputs.2d.Kolmogorov_ml ns.ml_ckpt_path=$ml_ckpt_path ns.inp_scale_factor=$inp_scale_factor ns.target_scale_factor=$target_scale_factor ns.expt_dir=$inference_dir ns.history=false ns.do_inference=1 ns.ml_correction=1 stop_interval=$stop_interval ns.ml_refinement_ratio=$ml_refinement_ratio amr.restart=$restart_ckpt amr.n_cell= 32 32 ns.cfl=0.4

#run baseline
./amr2d.gnu.haswell.MPI.ex inputs.2d.Kolmogorov_ml ns.ml_ckpt_path=$ml_ckpt_path ns.inp_scale_factor=$inp_scale_factor ns.target_scale_factor=$target_scale_factor ns.expt_dir=$inference_dir ns.history=false ns.do_inference=1 ns.ml_correction=0 stop_interval=$stop_interval ns.ml_refinement_ratio=$ml_refinement_ratio amr.restart=$restart_ckpt amr.n_cell= 32 32 ns.cfl=0.4



#build hdf5 file with results
h5_dir='/global/cscratch1/sd/jpathak/iamr_inference/kolmogorov/'$config_name$run_num
mkdir -p $h5_dir
module load python
source activate $JPATHAK/envs
#python analysis/cataloge.py
python analysis/inference_plots_to_h5.py $inference_dir $h5_dir 


