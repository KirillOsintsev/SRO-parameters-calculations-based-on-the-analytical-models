#!/bin/bash
# @Xin Liu, xin.liu@epfl.ch

set -o nounset # Treat unset variables as an error

/home/max/miniforge3/bin/python3 \
    calculate_1st_S-S_interaction_by_pressure_relaxation.py \
    /media/sf_host/SRO-parameters-calculator/data/input/potentials/average_potentials/Fe60Ni20Cr20.averaged.eam.alloy \
    X Fe Ni 5 3.52 eam/alloy > log12.txt

/home/max/miniforge3/bin/python3 \
    calculate_1st_S-S_interaction_by_pressure_relaxation.py \
    /media/sf_host/SRO-parameters-calculator/data/input/potentials/average_potentials/Fe60Ni20Cr20.averaged.eam.alloy \
    X Fe Fe 5 3.52 eam/alloy > log11.txt

/home/max/miniforge3/bin/python3 \
    calculate_1st_S-S_interaction_by_pressure_relaxation.py \
    /media/sf_host/SRO-parameters-calculator/data/input/potentials/average_potentials/Fe60Ni20Cr20.averaged.eam.alloy \
    X Ni Ni 5 3.52 eam/alloy > log22.txt
