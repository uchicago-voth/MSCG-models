#!/bin/bash

LAMMPS=${HOME}/bin/lmp_clx_UCG

SEED=${RANDOM}
NN=${SLURM_NTASKS}

ibrun -n 336 $LAMMPS -in input.init -var MID 1 -var SEED ${RANDOM} > ${SLURM_JOBID}.oe 2>&1
