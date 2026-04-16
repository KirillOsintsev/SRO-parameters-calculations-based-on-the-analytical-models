#!/bin/bash
# @Xin Liu, xin.liu@epfl.ch

set -o nounset # Treat unset variables as an error

python3 \
    calculate_2nd_S-S_interaction_by_pressure_relaxation.py \
    /media/sf_host/Project/data/input/potentials/average_potentials/Fe56Ni23Cr21.averaged.eam.alloy \
    X Cr Ni 5 3.52 eam/alloy > log12.txt

python3 \
    calculate_2nd_S-S_interaction_by_pressure_relaxation.py \
    /media/sf_host/Project/data/input/potentials/average_potentials/Fe56Ni23Cr21.averaged.eam.alloy \
    X Cr Cr 5 3.52 eam/alloy > log11.txt

python3 \
    calculate_2nd_S-S_interaction_by_pressure_relaxation.py \
    /media/sf_host/Project/data/input/potentials/average_potentials/Fe56Ni23Cr21.averaged.eam.alloy \
    X Ni Ni 5 3.52 eam/alloy > log22.txt
