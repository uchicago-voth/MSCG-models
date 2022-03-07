#!/bin/bash               

LAMMPS=${HOME}/bin/lmp_clx_UCG

SEED=${RANDOM}
#the molID for the Gag molecule that "controls" the cluster
#this is stochastic so chosen based on which Gag
#appears to be assembling after a few trial runs
MID=18

init=1

for i in {1..10}
do
    mkdir trial_${i}
    cd trial_${i}

    if [ $init -eq 1 ]
    then
	ln -s ../system.init
	ln -s ../bond.ff
	ln -s ../exclusion_withcoul.ff
	ln -s ../gaussoff.ff
	ln -s ../gausswall.ff
	ln -s ../cross_interactions.ff
	ln -s ../system.settings
	ln -s ../rates.txt
	ln -s ../contacts.txt
	ln -s ../input.ref
	ln -s ../eq0.data
    fi

    let OFFSET=($i-1)*336
    ibrun -n 336 -o $OFFSET $LAMMPS -in input.ref -var SEED ${RANDOM} > ${SLURM_JOBID}_${i}.oe 2>&1 &

    cd ..
done

wait

bash assess.bash
