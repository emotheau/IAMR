#!/bin/bash

for i in {8..20};do sbatch --dependency=singleton --job-name=generate_data_256 submit_job_256.sh $i; done
