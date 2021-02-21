#AMREX_HOME defines the directory in which we will find the BoxLib directory
AMREX_HOME ?= ../../../amrex

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION   = DOUBLE

DIM         = 2
COMP        = gnu

DEBUG       = TRUE
USE_MPI     = TRUE
USE_OMP     = FALSE
PROFILE     = FALSE

USE_CUDA = FALSE

USE_LIBTORCH = TRUE
LIBTORCH_DIR = /home/emotheau/GitCodes/libtorch_cpu

USE_SENSEI_INSITU = FALSE

EBASE = amr

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
