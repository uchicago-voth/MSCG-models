#!/bin/bash

path_to_perl_scripts=/Users/ajpak/Downloads/heteroENM

AVG=../../subset_CA.pdb
SUP=../../all_aligned_CA.pdb
MAP=../mapping.trans
NCG=35
NCUT=20.0
KBEND=0.5


perl ${path_to_perl_scripts}/fluc-match-str-pdb.pl $AVG $SUP $MAP bond.out fluct.out $NCUT $KBEND

perl ${path_to_perl_scripts}/fluc-match-8f.pl bond.out fluct.out $NCG 5000 0

# optional LAMMPS converter
perl ${path_to_perl_scripts}/convert_to_lammps.pl cgk.dat mass.dat cg.xyz henm.lammps.data
