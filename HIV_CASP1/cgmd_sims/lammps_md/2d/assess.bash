#!/bin/bash

LAMMPS=${HOME}/bin/lmp_clx_UCG
max_folder=""
max_cluster=0

SAVE_FOLDER=master_traj

# RENAME LAMMPS LOG
for i in {1..10}
do
    if [ -f trial_${i}/log.lammps ]
    then
	cp trial_${i}/log.lammps trial_${i}/log.1
    fi
done

# FIND TRIAL WITH LARGEST CLUSTER
fid=$( python assess_candidates.py )
max_cluster=$( grep "25000000  " trial_${fid}/log.1 | gawk '{print $21}' )
max_folder=trial_${fid}

echo "The candidate $max_folder has a final cluster size of $max_cluster"

# Now FIND NEXT SUFFIX INDEX FOR STORING
next_store_id=temp
if [ -d $SAVE_FOLDER ]
then
    cd $SAVE_FOLDER

    for i in {1..1000}
    do
	if [ ! -f eq${i}.lammpstrj ]
	then
	    next_store_id=${i}
	    break
	fi
    done

    cd ../
fi

echo "The candidate $max_folder will be saved with suffix id $next_store_id"

# SAVE FILES
if [ -d $max_folder ]
then
    cd $max_folder

    $LAMMPS -restart2data eq1.restart eq1.data -log log.convert
    cp eq1.data ../eq0.data
    
    cp log.1 ../${SAVE_FOLDER}/log.${next_store_id}
    cp dump1.lammpstrj ../${SAVE_FOLDER}/dump${next_store_id}.lammpstrj
    cp sub1.lammpstrj ../${SAVE_FOLDER}/sub${next_store_id}.lammpstrj
    cp eq1.lammpstrj ../${SAVE_FOLDER}/eq${next_store_id}.lammpstrj
    cp eq1.data ../${SAVE_FOLDER}/eq${next_store_id}.data

    cd ../
fi
