#!/bin/bash               

export OMP_NUM_THREADS=1

LAMMPS=${HOME}/bin/lmp_clx_UCG

SEED=${RANDOM}
NN=${SLURM_NTASKS}

init=1
#the molID for the Gag molecule that "controls" the cluster
#this is stochastic so chosen based on which Gag
#appears to be assembling after a few trial runs
MID=2321

for i in {1..10}
do
    mkdir trial_${i}
    cd trial_${i}

    if [ $init -eq 1 ]
    then
	ln -s ../system.init .
	ln -s ../lipid.ff .
	ln -s ../bond.ff .
	ln -s ../exclusion.ff .
	ln -s ../gaussoff.ff .
	ln -s ../gausswall.ff .
	ln -s ../cross_interactions.ff .
	ln -s ../system.settings .
	ln -s ../eq0.data .
	ln -s ../rates.txt .
	ln -s ../contacts.txt .
	ln -s ../input.run .
    fi

    let OFFSET=($i-1)*784
    ibrun -n 784 -o $OFFSET $LAMMPS -in input.run -var SEED ${RANDOM} -var MID ${MID} > ${SLURM_JOBID}_${i}.oe 2>&1 &

    cd ../
done

wait

bash assess.bash
