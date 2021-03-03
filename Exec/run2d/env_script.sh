#!/bin/bash

export AMREX_HOME=/project/projectdirs/dasrepo/jpathak/amrex
module load pytorch
module switch PrgEnv-intel PrgEnv-gnu
module load cray-mpich
module load cray-mpich-abi
module load cmake


