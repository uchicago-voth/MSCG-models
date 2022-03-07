#!/bin/bash

i=35

EDCG=/home/ajpak/CGMGC/build/cgmgc #serial

STDOUT='stdout.log'
STDERR='stderr.log'

# 224 CA/residue per chain
ncg=${i}
input_psf="../subset.psf"
input_dcd="../all_aligned.dcd"
controlLog="controlLog.out"
fullLog="controlLog.out"
cgtraj="minCgtraj.xyz"
lambda=0
step=10

$EDCG --linear-model \
    --num-CG-sites $ncg \
    --input-psf $input_psf \
    --input-dcd $input_dcd \
    --model-log    $controlLog    \
    --num-cov-dof   -1            \
    --num-models  1           \
    --num-res-bs-reps 0 \
    --num-res-bs-models 0 \
    --num-bs-reps 0        \
    --sim-ann-acceptance 0.97        \
    --sim-ann-control    0.07       \
    --use-bca-ci 0        \
    --stride $step   >$STDOUT 2>$STDERR

