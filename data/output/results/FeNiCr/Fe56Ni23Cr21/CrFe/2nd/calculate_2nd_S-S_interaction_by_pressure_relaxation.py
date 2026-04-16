#!/usr/bin/env python3
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
    template.template += "region sphere sphere 0 0 0 0.05 units box\n"
    template.template += "group solute region sphere\n"
    template.template += "set group solute type 2\n"

    template.template += "region sphere2 sphere 1 0 0 0.05 units lattice\n"
    template.template += "group solute2 region sphere2\n"
    template.template += "set group solute2 type 3\n"

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
            safe_template_string += line + "\n"
    return CustomTemplate(safe_template_string)

class CustomTemplate(Template):
    delimiter = "?"

lammps_input_template = safe_template_string("""
units           metal
boundary        p p p
lattice         fcc ?lattice_parameter
region          simbox block 0 ?ncells 0 ?ncells 0 ?ncells units lattice
create_box      3 simbox
create_atoms    1 box
pair_style      eam/alloy
pair_coeff      * * ?potential ?matrix_type ?solute_type ?solute_type2

variable        p equal press
variable        pe equal pe
thermo 25
thermo_style    custom step temp pe press fnorm
""")

if __name__ == "__main__":
    main()
