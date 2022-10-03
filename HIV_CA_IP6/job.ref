#!/bin/bash               

#SBATCH -J hiv_r15
#SBATCH -o hiv1_l.out
#SBATCH -e hiv1_l.err
#SBATCH -p development
#SBATCH -N 2
#SBATCH -n 112
#SBATCH -t 00:20:00
#SBATCH -A CHE20010

export OMP_NUM_THREADS=1

module unload intel
module unload impi
module unload mvapich2-x
module load intel/19.0.5
module load impi
module load fftw3
module load gsl

LAMMPS=/work2/07732/tg870312/frontera/lammps-21Jul20/src/lmp_intel_cpu_intelmpi

ibrun -n 112 $LAMMPS -in input -var SEED $RANDOM

