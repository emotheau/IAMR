#!/bin/bash

module load parallel

seq 0 8 | srun parallel --jobs 8 write_chunked_h5.py {}
