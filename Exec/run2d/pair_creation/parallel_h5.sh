#!/bin/bash

module load parallel

res=256
user_ratio=4
processed_files_dir='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'$res'/ur'$user_ratio

for idx in 0 3 4 5 6 7 8 9 10 11 12 13 14
do
  num_files=$(find $processed_files_dir/$idx/plots/* -maxdepth 0 -type d | wc -l)
  echo $num_files
  num_chunks=$(($num_files / 1000))
  seq 0 $(($num_chunks-1)) | srun parallel --jobs $num_chunks python hdf5_writer.py {} $idx
#  seq 0 $(($num_chunks-1)) | srun parallel --jobs $num_chunks echo {} $idx
  sleep 10
done


#seq 0 6 | srun parallel --jobs 7 python write_chunked_with_history_h5.py {} 3
#
#sleep 10
#
#seq 0 7 | srun parallel --jobs 8 python write_chunked_with_history_h5.py {} 4
#
#sleep 10
#
#seq 0 6 | srun parallel --jobs 7 python write_chunked_with_history_h5.py {} 5
#
#sleep 10
#
#seq 0 25 | srun parallel --jobs 26 python write_chunked_with_history_h5.py {} 0
