#!/bin/bash

module unload intel
module unload impi
module load intel/18.0.5
module load mvapich2
module load fftw3
module load gsl
module load python2

LAMMPS=~/bin/lmp_clx_UCG
NEWREM=~/bin/rangefinder_no_gro.x

for i in {1..500}
do
    let j=${i}-1
    OLD=iter${j}
    NEW=iter${i}

    mkdir $NEW
    mkdir ${NEW}/RF

    cd $NEW

    ln -s ../bond.ff ./
    ln -s ../exclusion.ff ./
    ln -s ../system.init ./
    ln -s ../system.settings ./
    ln -s ../system.data ./
    ln -s ../input.rem ./

    cp ../${OLD}/RF/gausswall.ff ./

    mpirun -np 16 $LAMMPS -in input.rem -var SEED ${RANDOM}

    if [ ! -e ./eq1.data ]
    then
	exit 1
    fi
    
    cd RF
    ln -s ../../init/RF/top.in ./
    ln -s ../../init/RF/control.in ./

    #ulimit -S -n 2048
    $NEWREM -l ../dump1.lammpstrj    
    if [ ! -e ./25_25.hist ]
    then
	exit 1
    fi

    # learning rate schedule
    lr=1.0
    if [ $i -lt 20 ]
    then
	lr=0.80
    if [ $i -lt 50 ]
    then
	lr=0.20
    elif [ $i -lt 100 ]
    then
	lr=0.10
    elif [ $i -lt 150 ]
    then
	lr=0.01
    elif [ $i -lt 200 ]
    then
	lr=0.001
    elif [ $i -lt 500 ]
    then
	lr=0.0001
    fi

    python ../../evaluate_and_update_gauss_rem.py $lr

    if [ ! -e ./gausswall.ff ]
    then
	exit 1
    fi

    rm *.dist *.hist
    rm ../dump1.lammpstrj

    cd ../../

done

