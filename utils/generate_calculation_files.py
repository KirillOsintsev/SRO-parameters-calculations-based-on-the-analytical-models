#!/usr/bin/env python3
import os
import json
from pathlib import Path

def create_shell_script(pair_path, nn_type, config, pair):
    """Create the run_calculation shell script for a specific pair and nearest neighbor."""
    script_template = '''#!/bin/bash
# @Xin Liu, xin.liu@epfl.ch

set -o nounset # Treat unset variables as an error

python3 \\
    calculate_{nn}_S-S_interaction_by_pressure_relaxation.py \\
    /media/sf_host/Project/data/input/potentials/average_potentials/{alloy_name}.averaged.eam.alloy \\
    X {element1} {element2} {ncells} {lattice_constant} {potential_style} > log12.txt

python3 \\
    calculate_{nn}_S-S_interaction_by_pressure_relaxation.py \\
    /media/sf_host/Project/data/input/potentials/average_potentials/{alloy_name}.averaged.eam.alloy \\
    X {element1} {element1} {ncells} {lattice_constant} {potential_style} > log11.txt

python3 \\
    calculate_{nn}_S-S_interaction_by_pressure_relaxation.py \\
    /media/sf_host/Project/data/input/potentials/average_potentials/{alloy_name}.averaged.eam.alloy \\
    X {element2} {element2} {ncells} {lattice_constant} {potential_style} > log22.txt
'''

    script_content = script_template.format(
        nn=nn_type,
        alloy_name=config['alloy_name'],
        element1=pair[0],
        element2=pair[1],
        ncells=config['ncells'],
        lattice_constant=config['lattice_constant'],
        potential_style=config['potential_style']
    )
    
    script_path = pair_path / nn_type / f"run_calculation_{nn_type}.sh"
    
    with open(script_path, 'w') as f:
        f.write(script_content)
        
    return script_path
    
    # Make the script executable
    os.chmod(script_path, 0o755)

def create_python_script(pair_path, nn_type, config):
    """Create the Python calculation script for a specific pair and nearest neighbor."""
    # Define nearest neighbor coordinates based on crystal structure
    nn_coords = {
        'bcc': {
            '1st': '0.5 0.5 0.5',      # First nearest neighbor for BCC
            '2nd': '0 0 1',            # Second nearest neighbor for BCC
            '3rd': '0 1 1',            # Third nearest neighbor for BCC
            '4th': '0.5 0.5 1.5'       # Fourth nearest neighbor for BCC
        },
        'fcc': {
            '1st': '0.5 0.5 0',        # First nearest neighbor for FCC
            '2nd': '1 0 0',            # Second nearest neighbor for FCC
            '3rd': '1 1 0',            # Third nearest neighbor for FCC
            '4th': '1.5 0.5 0'         # Fourth nearest neighbor for FCC
        }
    }

    template = '''#!/usr/bin/env python3
"""Calculate minimum energy configuration using pressure relaxation.

Parameters
----------
potential : str
    Path to the EAM potential table (Lammps setfl format)
matrix_type : str
solute_type : str
ncells : str
    number of unit cells along the directions
lattice_parameter : str
    guess for the lattice parameter
"""
import sys
from string import Template
from scipy import optimize
from lammps import lammps
from mpi4py import MPI

def main():
    potential = sys.argv[1]
    matrix_type = sys.argv[2]
    solute_type = sys.argv[3]
    solute_type2 = sys.argv[4]
    ncells = sys.argv[5]
    lattice_parameter = sys.argv[6]
    template = safe_template_string(
        lammps_input_template.safe_substitute(
            matrix_type=matrix_type,
            solute_type=solute_type,
            solute_type2=solute_type2,
            potential=potential,
            ncells=ncells
        )
    )

    # Find lattice parameter at zero pressure in configuration with solute
    if MPI.COMM_WORLD.rank == 0:
        print("Find lattice parameter at zero pressure in configuration with solute")
    template.template += "region sphere sphere 0 0 0 0.05 units box\\n"
    template.template += "group solute region sphere\\n"
    template.template += "set group solute type 2\\n"

    template.template += "region sphere2 sphere {sphere2_coords} 0.05 units lattice\\n"
    template.template += "group solute2 region sphere2\\n"
    template.template += "set group solute2 type 3\\n"

    a = float(lattice_parameter) * 0.95
    b = float(lattice_parameter) * 1.1
    x1, r = optimize.brentq(
        objective_function, a, b, args=(template, False),
        full_output=True
    )
    print("#########")
    residual_pressure = objective_function(x1, template, file="data.1")
    if MPI.COMM_WORLD.rank == 0:
        report(x1, residual_pressure, r)

def report(x, y, r):
    print("lattice parameter (Angstroem): {:.15f}".format(x))
    print("residual pressure (bar):      {:.15f}".format(y))
    print("convergence info:")
    print(r)

def objective_function(lattice_parameter, template, verbose=False, file=None):
    x = "{:.15f}".format(lattice_parameter)
    lammps_input = template.substitute(lattice_parameter=x)
    if verbose:
        cmdargs = []
    else:
        cmdargs = ["-echo", "none", "-screen", "none"]
    lmp = lammps(comm=MPI.COMM_WORLD, cmdargs=cmdargs)
    for line in lammps_input.splitlines():
        lmp.command(line)
    lmp.command("min_style cg")
    lmp.command("minimize 0.0 1.0e-6 100000 100000")
    lmp.command("min_style fire")
    lmp.command("minimize 0.0 1.0e-6 100000 100000")
    pressure = lmp.extract_variable('p', 'all', 0)
    energy = lmp.extract_variable('pe', 'all', 0)
    if MPI.COMM_WORLD.rank == 0:
        print(pressure)
        print(energy)
    if file is not None:
        lmp.command("write_data {:s}".format(file))
    lmp.close()
    return pressure

def safe_template_string(multiline_string):
    safe_template_string = ""
    for line in multiline_string.splitlines():
        if line.rstrip().endswith("&"):
            safe_template_string += line.rstrip().rstrip("&")
        else:
            safe_template_string += line + "\\n"
    return CustomTemplate(safe_template_string)

class CustomTemplate(Template):
    delimiter = "?"

lammps_input_template = safe_template_string("""
units           metal
boundary        p p p
lattice         {crystal_type} ?lattice_parameter
region          simbox block 0 ?ncells 0 ?ncells 0 ?ncells units lattice
create_box      3 simbox
create_atoms    1 box
pair_style      {potential_style}
pair_coeff      * * ?potential ?matrix_type ?solute_type ?solute_type2

variable        p equal press
variable        pe equal pe
thermo 25
thermo_style    custom step temp pe press fnorm
""")

if __name__ == "__main__":
    main()
'''

    # Replace placeholders in the template
    script_content = template.replace('{crystal_type}', config['phase'])
    script_content = script_content.replace('{potential_style}', config['potential_style'])
    script_content = script_content.replace('{sphere2_coords}', nn_coords[config['phase']][nn_type])
    
    script_path = pair_path / nn_type / f"calculate_{nn_type}_S-S_interaction_by_pressure_relaxation.py"
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make the script executable
    os.chmod(script_path, 0o755)

def main():
    # Load configuration
    try:
        config_files = list(Path('.').glob('*/{alloy_name}.config.json'))
        if not config_files:
            raise FileNotFoundError("No config.json found in subdirectories")
        
        config_path = config_files[0]
        with open(config_path, 'r') as f:
            config = json.load(f)
    except Exception as e:
        print(f"Error loading configuration: {e}")
        return

    alloy_path = Path(config['alloy_name'])
    
    # Generate calculation files for each pair and nearest neighbor
    for pair in config['pairs']:
        pair_path = alloy_path / ''.join(pair)
        for nn in config['nearest_neighbors']:
            print(f"Generating calculation files for {pair[0]}-{pair[1]} {nn}")
            create_shell_script(pair_path, nn, config, pair)
            create_python_script(pair_path, nn, config)

if __name__ == "__main__":
    main()