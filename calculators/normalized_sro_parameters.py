#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:35:24 2024

@author: osintsevkirill

References
----------
The normalized SRO metrics computed here (β_i, β_0, and the reported FCC/BCC averages)
follow formulas used in the solute-strengthening with SRO framework:

Nag, S. and Curtin, W. A. (2024).
"Solute-strengthening in metal alloys with short-range order", Acta Materialia 263, 119472.
DOI: 10.1016/j.actamat.2023.119472
"""

import os
import re


def parse_misfit_volume_file(file_name):
    """Parse the misfit volume results file."""
    with open(file_name, 'r') as file:
        lines = file.readlines()
    
    # Extract misfit volumes
    misfit_volumes = {}
    misfit_index = lines.index("Element\tMisfit Volume (Angstroem^3)\n") + 1
    for line in lines[misfit_index:]:
        if line.strip():
            element, volume = line.split()
            misfit_volumes[element.strip()] = float(volume.strip())
    
    return misfit_volumes


def calculate_normalized_sro(phase, composition, misfit_volumes, sro_parameters):
    """Calculate normalized SRO parameters for all shells and compute averages for FCC/BCC phases."""
    # Ensure compositions are normalized to fractions
    total_composition = sum(composition.values())
    normalized_composition = {
        element: comp / total_composition for element, comp in composition.items()
    }

    # Calculate beta_0
    beta_0 = sum(
        comp * misfit_volumes[element]**2
        for element, comp in normalized_composition.items()
    )

    # Store results
    beta_values = {}
    averages = {}

    for temperature, shells in sro_parameters.items():
        beta_shells = {}

        # Calculate beta_i for each shell
        for shell_idx in ["1st", "2nd", "3rd", "4th"]:
            if shell_idx not in shells:
                continue  # Skip if the shell is not present in the data

            alpha_values = shells[shell_idx]
            beta_i = 0

            # Iterate through element pairs
            elements = list(normalized_composition.keys())
            for i, elem_n in enumerate(elements):
                for elem_m in elements[i + 1:]:
                    delta_v_diff = misfit_volumes[elem_n] - misfit_volumes[elem_m]
                    alpha_key = f"alpha_{elem_n}{elem_m}" if f"alpha_{elem_n}{elem_m}" in alpha_values else f"alpha_{elem_m}{elem_n}"

                    if alpha_key in alpha_values:
                        beta_i += (
                            normalized_composition[elem_n] * normalized_composition[elem_m] *
                            (delta_v_diff ** 2) * alpha_values[alpha_key]
                        )

            beta_shells[shell_idx] = beta_i
            #print(beta_shells)

        beta_values[temperature] = beta_shells
        #print(beta_values)

        # Calculate averages based on the phase
        if phase == "fcc":
            beta_1_beta_0 = beta_shells.get("1st", 0) / beta_0
            beta_2 = beta_shells.get("2nd", 0)
            beta_3 = beta_shells.get("3rd", 0)
            beta_4 = beta_shells.get("4th", 0)
            averages[temperature] = {
                "beta_average/beta_0": (beta_2 + 4 * beta_3 + 2 * beta_4) / beta_0,
                "beta_1/beta_0": beta_1_beta_0
            }
            

        elif phase == "bcc":
            beta_1 = beta_shells.get("1st", 0)
            beta_2 = beta_shells.get("2nd", 0)
            beta_3 = beta_shells.get("3rd", 0)
            beta_4 = beta_shells.get("4th", 0)
            averages[temperature] = {
                "beta_1_average/beta_0": (1.5 * beta_1 + beta_2) / beta_0,
                "beta_2_average/beta_0": (beta_3 + 1.8 * beta_4) / beta_0,
            }

    return beta_0, beta_values, averages


def write_output(alloy_path, alloy_name, phase, composition, beta_0, beta_values, averages):
    """Write results to output file."""
   
    output_file = alloy_path / f"{alloy_name}_normalized_sro_parameters.txt"
    
    
    with open(output_file, 'w') as file:
        # Write composition
        file.write(f"Alloy Composition:\n")
        for element, comp in composition.items():
            file.write(f"{element}: {comp:.2f}\n")
        file.write("\n")
        
        # Write phase and beta_0
        file.write(f"Phase: {phase.upper()}\n\n")
        file.write(f"β_0: {beta_0:.6f}\n\n")
        
        # Write SRO parameters for each temperature
        file.write("Temperature\tDistance\tβ_i\tNormalized SRO (β_i/β_0)\n")
        for temperature, distances in beta_values.items():
            for distance, beta_i in distances.items():
                normalized_sro = beta_i / beta_0  # Calculate normalized SRO
                file.write(f"{temperature}K\t{distance}\t{beta_i:.6f}\t{normalized_sro:.6f}\n")
        
        file.write("\n")
        
        # Write averages for each temperature
        file.write("Temperature\tAverage SRO Metrics\n")
        for temperature, avg_values in averages.items():
            file.write(f"{temperature}K:\n")
            for key, value in avg_values.items():
                file.write(f"  {key}: {value:.6f}\n")
        # Also write average β_1/β_0 for FCC phase
            if phase == "fcc":
                beta_1_b0_avg = beta_values.get(temperature, {}).get("1st", 0) / beta_0  # Calculate β_1/β_0 average
                file.write(f"  β_1/β_0: {beta_1_b0_avg:.17f}\n")        
    
    print(f"Results written to {output_file}")



def calculate_normalized_sro_parameters(sro_parameters, sro_output_file, misfit_volume_file, phase, composition, alloy_path, alloy_name):
    
    # Parse files
    misfit_volumes = parse_misfit_volume_file(misfit_volume_file)
    
    # Calculate normalized SRO parameters
    beta_0, results, averages = calculate_normalized_sro(phase, composition, misfit_volumes, sro_parameters)
    
    # Write results to output file
    write_output(alloy_path, alloy_name, phase, composition, beta_0, results, averages)
    
    return averages
  
   
   
    
     
