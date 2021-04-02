#!/bin/bash

for i in {1..20};do sbatch --dependency=singleton --job-name=generate_data submit_job.sh $i; done
