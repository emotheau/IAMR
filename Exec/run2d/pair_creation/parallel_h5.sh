#!/bin/bash

module load parallel

seq 0 6 | srun parallel --jobs 7 python write_chunked_with_history_h5.py {} 3

sleep 10

seq 0 7 | srun parallel --jobs 8 python write_chunked_with_history_h5.py {} 4

sleep 10

seq 0 6 | srun parallel --jobs 7 python write_chunked_with_history_h5.py {} 5

sleep 10

seq 0 25 | srun parallel --jobs 26 python write_chunked_with_history_h5.py {} 0
