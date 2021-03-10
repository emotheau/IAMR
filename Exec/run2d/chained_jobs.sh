#!/bin/bash

#jobid1=$(sbatch --parsable submit_job.sh)
#jobid2=$(sbatch --parsable --dependency=after:$jobid1 submit_job.sh)
#jobid3=$(sbatch --parsable --dependency=after:$jobid2 submit_job.sh)
#sbatch --dependency=after:$jobid3 submit_job.sh

for i in {0..5};do sbatch --dependency=singleton --job-name=generate_data submit_job.sh ; done
