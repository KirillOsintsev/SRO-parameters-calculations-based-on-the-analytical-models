#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes / References
------------------
This module computes misfit volumes relative to a composition-weighted average atomic volume.
The linear mixture rule is consistent with (a Vegard-like) approximation commonly used in
solute-strengthening models, e.g.:

Varvenne, C., Luque, A., Curtin, W. A. (2016).
"Theory of strengthening in fcc high entropy alloys", Acta Materialia 118, 164–176.
DOI: 10.1016/j.actamat.2016.07.040
"""
import json
import re

def load_database(db_path):
    """
    Load the alloy and element data from the database file.
    :param db_path: Path to the JSON database file.
    :return: Parsed database as a dictionary.
    """
    with open(db_path, 'r') as file:
        return json.load(file)


def extract_atomic_volumes(database, elements):
    """
    Extract atomic volumes for specific elements from the database.
    :param database: Dictionary containing element data.
    :param elements: List of elements to extract atomic volumes for.
    :return: Dictionary of atomic volumes.
    """
    atomic_volumes = {}
    for element in elements:
        if element in database["elements"]:
            atomic_volumes[element] = database["elements"][element]["Vn"]
        else:
            raise ValueError(f"Element {element} not found in the database.")
    return atomic_volumes


def calculate_misfit_volumes(atomic_volumes, concentrations):
    """
    Calculate misfit volumes for all elements in the alloy.
    :param atomic_volumes: Dictionary of atomic volumes for all elements.
    :param concentrations: List of concentrations for each element.
    :return: Dictionary of misfit volumes.
    """
    average_volume = sum(c * atomic_volumes[element] for c, element in zip(concentrations, atomic_volumes.keys()))
    misfit_volumes = {}

    for element, volume in atomic_volumes.items():
        misfit_volumes[element] = volume - average_volume

    return average_volume, misfit_volumes


def parse_composition(composition_str):
    """
    Parse the alloy composition string into elements and their concentrations.
    :param composition_str: Composition string in the format "Nb2Co32Cr24Ni40Mn2".
    :return: Tuple of element list and concentration list.
    """
    pattern = r"([A-Z][a-z]*)(\d*)"
    matches = re.findall(pattern, composition_str)

    elements = []
    concentrations = []
    total_atoms = 0

    for element, count in matches:
        count = int(count) if count else 1
        elements.append(element)
        concentrations.append(count)
        total_atoms += count

    # Normalize concentrations to percentages
    concentrations = [c / total_atoms for c in concentrations]

    return elements, concentrations


def save_results(misfit_volumes, concentrations, output_file, average_volume):
    """
    Save misfit volumes and alloy composition to an output file.
    :param misfit_volumes: Dictionary of misfit volumes.
    :param concentrations: List of concentrations for each element.
    :param output_file: Path to the output file.
    """
    with open(output_file, 'w') as file:
        file.write("Alloy Composition:\n")
        for element, concentration in zip(misfit_volumes.keys(), concentrations):
            file.write(f"{element}: {concentration:.2f}\n")
        file.write(f"Average volume: {average_volume:.3f}\n")
        file.write("\nElement\tMisfit Volume (Angstroem^3)\n")
        for element, misfit in misfit_volumes.items():
            file.write(f"{element}\t{misfit:.6f}\n")  # Fixed format with six decimal places
    print(f"Misfit volumes and composition saved to {output_file}")


def misfit_volume_calculator(database_path, alloy_path, alloy_name, composition, components, concentrations):
    db_path = database_path  # Path to the JSON database file

    # Generate output file name based on composition
    output_file = alloy_path / f"{alloy_name}_misfit_volumes_results.txt"

    # Load database
    database = load_database(db_path)

    # Extract atomic volumes for specified elements
    atomic_volumes = extract_atomic_volumes(database, components)

    # Calculate misfit volumes
    average_volume, misfit_volumes = calculate_misfit_volumes(atomic_volumes, concentrations)
    

    # Save results including composition
    save_results(misfit_volumes, concentrations, output_file, average_volume)


