#!/usr/bin/env python3
#
# Third-party attribution
# -----------------------
# This file contains code adapted from `average_atom.py` from the `matscipy` project:
# - Copyright 2021 Lars Pastewka (U. Freiburg)
# - Copyright 2020 Wolfram G. Nöhring (U. Freiburg)
# Source: https://github.com/libAtoms/matscipy
#
# The original code is distributed under the GNU General Public License (GPL),
# version 2 or (at your option) any later version. See <http://www.gnu.org/licenses/>.
import os
import json
from collections import namedtuple
from pathlib import Path
import numpy as np
from .composition import get_composition_from_user, get_alloy_name_and_concentrations


"""Conversion factor from Hartree to Electronvolt. \
Use 1 Hartree = 27.2 eV to be consistent with Lammps EAM."""
_hartree_in_electronvolt = 27.2

"""Conversion factor from Bohr radii to Angstrom. \
Use 1 Bohr radius = 0.529 A to be consistent with Lammps EAM."""
_bohr_in_angstrom = 0.529

"""Conversion factor for charge function of EAM potentials \
of kind 'eam' (DYNAMO funcfl format)."""
_hartree_bohr_in_electronvolt_angstrom = _hartree_in_electronvolt * _bohr_in_angstrom

# Todo: replace by data class (requires Python > 3.7)
class EAMParameters(
    namedtuple(
        "EAMParameters",
        "symbols atomic_numbers "
        "atomic_masses lattice_constants crystal_structures "
        "number_of_density_grid_points "
        "number_of_distance_grid_points "
        "density_grid_spacing distance_grid_spacing "
        "cutoff",
    )
):
    """Embedded Atom Method potential parameters

    :param array_like symbols: Symbols of the elements coverered by
        this potential (only for eam/alloy and 
        eam/fs, EMPTY for eam
    :param array_like atomic_numbers: Atomic numbers of the elements 
        covered by this potential
    :param array_like atomic_masses: Atomic masses of the elements 
        covered by this potential
    :param array_like lattice_constants: Lattice constant of a pure crystal 
        with crystal structure as specified in crystal_structures
    :param array_like crystal_structures: Crystal structure of the pure metal.
    :param int number_of_density_grid_points: Number of grid points 
        of the embedding energy functional
    :param int number_of_distance_grid_points: Number of grid points of 
        the electron density function and the pair potential
    :param float density_grid_spacing: Grid spacing in electron density space
    :param float distance_grid_spacing: Grid spacing in pair distance space
    :param float cutoff: Cutoff distance of the potential
    """

    __slots__ = ()

###

def _strip_comments_from_line(string, marker="#"):
    """Strip comments from lines but retain newlines

    Parameters
    ----------
    string : str
        string which may contain comments
    marker : str
        marker which indicates the start of a comment

    Returns
    -------
    stripped_string : str
        string without commments; if the string terminated
        with a newline, then the newline is retained
    """
    start = string.find(marker)
    if start != -1:
        stripped_string = string[:start]
        if string.endswith("\n"):
            stripped_string += "\n"
    else:
        stripped_string = string
    return stripped_string

def read_eam(eam_file, kind="eam/alloy"):
    """Read a tabulated EAM potential
    
    There are differnt flavors of EAM, with different storage
    formats. This function supports a subset of the formats supported
    by Lammps (http://lammps.sandia.gov/doc/pair_eam.html),
    * eam (DYNAMO funcfl format)
    * eam/alloy (DYNAMO setfl format)
    * eam/fs (DYNAMO setfl format)
    
    Parameters
    ----------
    eam_file : string
        eam alloy file name 
    kind : {'eam', 'eam/alloy', 'eam/fs'}
        kind of EAM file to read

    Returns
    -------
    source : string
        Source informations or comment line for the file header
    parameters : EAMParameters
        EAM potential parameters
    F : array_like
        contain the tabulated values of the embedded functions
        shape = (nb elements, nb of data points)
    f : array_like
        contain the tabulated values of the density functions
        shape = (nb elements, nb of data points)
    rep : array_like
        contain the tabulated values of pair potential
        shape = (nb elements,nb elements, nb of data points)
    """
    supported_kinds = ["eam", "eam/alloy", "eam/fs"]
    if kind not in supported_kinds: 
        raise ValueError(f"EAM kind {kind} not supported")
    # with open(eam_file, 'r', encoding='utf-8') as file:
    #     eam = file.readlines()

    if kind == "eam":
        with open(eam_file, 'r', encoding='utf-8') as file:
            # ignore comment characters on first line but strip them from subsequent lines
            lines = [file.readline()]
            lines.extend(_strip_comments_from_line(line) for line in file.readlines())
        # reading first comment line as source for eam potential data
        source = lines[0].strip()

        words = lines[1].strip().split()
        if len(words) != 4:
            raise ValueError(
                "expected four values on second line of EAM setfl file: "
                "atomic number, mass, lattice constant, lattice"
            )
        # Turn atomic numbers and masses, lattice parameter, and crystal
        # structures into arrays to be consistent with the other EAM styles
        atomic_numbers = np.array((int(words[0]), ), dtype=int)
        atomic_masses = np.array((float(words[1]), ), dtype=float)
        lattice_parameters = np.array((float(words[2]),), dtype=float)
        crystal_structures = np.empty(1).astype(str)
        crystal_structures[0] = words[3] 

        words = lines[2].strip().split()
        if len(words) != 5:
            raise ValueError(
                "expected five values on third line of EAM setfl file: "
                "Nrho, drho, Nr, dr, cutoff"
            )
        Nrho = int(words[0])       # Nrho (number of values for the embedding function F(rho))
        drho = float(words[1])     # spacing in density space
        Nr = int(words[2])         # Nr (number of values for the effective charge function Z(r) and density function rho(r))
        dr = float(words[3])       # spacing in distance space
        cutoff = float(words[4])
        parameters = EAMParameters(
            np.zeros(1), atomic_numbers, atomic_masses, lattice_parameters, 
            crystal_structures, Nrho, Nr, drho, dr, cutoff
        )

        # Strip empty lines 
        remaining_lines = [line for line in lines[3:] if len(line.strip()) > 0]
        remaining_words = []
        for line in remaining_lines:
            words = line.split()
            remaining_words.extend(words)
        expected_length = Nrho + 2 * Nr
        true_length = len(remaining_words)
        if true_length != expected_length:
            raise ValueError(f"expected {expected_length} tabulated values, but there are {true_length}")
        data = np.array(remaining_words, dtype=float)
        F = data[0:Nrho]
        # 'eam' (DYNAMO funcfl) tables contain the charge function :math:`Z`,
        # and not the pair potential :math:`\phi`. :math:`Z` needs to be 
        # converted into :math:`\phi` first, which involves unit conversion.
        # To be consistent with the other eam styles (and avoid complications
        # later), we convert into :math:`r*\phi`, where :math:`r` is the pair distance, i.e.
        # r = np.arange(0, rep.size) * dr
        charge = data[Nrho:Nrho+Nr]
        rep = charge**2
        rep *= _hartree_bohr_in_electronvolt_angstrom 
        f = data[Nrho+Nr:2*Nr+Nrho]
        # Reshape in order to be consistent with other EAM styles
        return source, parameters, F.reshape(1, Nrho), f.reshape(1, Nr), rep.reshape(1, 1, Nr)

    if kind in ["eam/alloy", "eam/fs"]:
        """eam/alloy and eam/fs have almost the same structure, except for the electron density section"""
        with open(eam_file, 'r', encoding='utf-8') as file:
            # ignore comment characters on first line but strip them from subsequent lines
            lines = [file.readline() for _ in range(3)]
            lines.extend(_strip_comments_from_line(line) for line in file.readlines())
        # reading 3 first comment lines as source for eam potential data
        source = "".join(line.strip() for line in lines[:3])

        words = lines[3].strip().split()
        alleged_num_elements = int(words[0])
        elements = words[1:]
        true_num_elements = len(elements)
        if alleged_num_elements != true_num_elements:
            raise ValueError(
                f"Header claims there are tables for {alleged_num_elements} elements, "
                f"but actual element list has {true_num_elements} elements: {' '.join(elements)}"
            )

        words = lines[4].strip().split()
        Nrho = int(words[0])     # Nrho (number of values for the embedding function F(rho))
        drho = float(words[1])   # spacing in density space
        Nr = int(words[2])       # Nr (number of values for the effective charge function Z(r) and density function rho(r))
        dr = float(words[3])     # spacing in distance space
        cutoff = float(words[4])

        # Strip empty lines and check that the table contains the expected number of values
        remaining_lines = [line for line in lines[5:] if len(line.strip()) > 0]
        remaining_words = []
        for line in remaining_lines:
            words = line.split()
            remaining_words.extend(words)
        if kind == "eam/fs":
            expected_num_density_functions_per_element = true_num_elements
        else:
            expected_num_density_functions_per_element = 1
        expected_num_words_per_element = (
            4 + 
            Nrho + 
            expected_num_density_functions_per_element * Nr 
        )
        expected_num_pair_functions = np.sum(np.arange(1, true_num_elements+1)).astype(int)
        expected_length = true_num_elements * expected_num_words_per_element + expected_num_pair_functions * Nr
        true_length = len(remaining_words)
        if true_length != expected_length:
            raise ValueError(f"expected {expected_length} tabulated values, but there are {true_length}")

        atomic_numbers = np.zeros(true_num_elements, dtype=int)
        atomic_masses = np.zeros(true_num_elements)
        lattice_parameters = np.zeros(true_num_elements)
        crystal_structures = np.empty(true_num_elements).astype(str) # fixme: be careful with string length
        F = np.zeros((true_num_elements, Nrho))
        for i in range(true_num_elements):
            offset = i * expected_num_words_per_element
            atomic_numbers[i] = int(remaining_words[offset])
            atomic_masses[i] = float(remaining_words[offset+1])
            lattice_parameters[i] = float(remaining_words[offset+2])
            crystal_structures[i] = remaining_words[offset+3]
            F[i, :] = np.array(remaining_words[offset+4:offset+4+Nrho], dtype=float)

        # Read data for individual elemements
        if kind == "eam/alloy":
            f = np.zeros((true_num_elements, Nr))
            for i in range(true_num_elements):
                offset = i * expected_num_words_per_element + 4 + Nrho
                f[i, :] = np.array(remaining_words[offset:offset+Nr], dtype=float)
        if kind == "eam/fs":
            f = np.zeros((true_num_elements, true_num_elements, Nr))
            for i in range(true_num_elements):
                offset = i * expected_num_words_per_element + 4 + Nrho
                for j in range(true_num_elements):
                    f[i, j, :] = np.array(remaining_words[offset+j*Nr:offset+(j+1)*Nr], dtype=float)

        # Read pair data
        rep = np.zeros((true_num_elements, true_num_elements, Nr))
        rows, cols = np.tril_indices(true_num_elements)
        for pair_number, (i, j) in enumerate(zip(rows, cols)):
            offset = true_num_elements * expected_num_words_per_element + pair_number * Nr
            rep[i, j, :] = np.array(remaining_words[offset:offset+Nr], dtype=float)
            rep[j, i, :] = rep[i, j, :]

        parameters = EAMParameters(
            elements, atomic_numbers, atomic_masses, 
            lattice_parameters, crystal_structures, 
            Nrho, Nr, drho, dr, cutoff
        )
        return source, parameters, F, f, rep

def write_eam(source, parameters, F, f, rep, out_file, kind="eam"):
    """Write an eam lammps format file
    # ... (docstring) ...
    """
    # Unpack parameters for the potential being written (includes N+1 elements now)
    elements, atomic_numbers, atomic_masses, lattice_parameters, crystal_structures = parameters[:5]
    Nrho, Nr, drho, dr, cutoff = parameters[5:]
    num_elements = len(elements) # This will be N+1 for averaged potentials

    if kind == "eam":
         # --- This block needs fixing if used ---
        print(f"Warning: Writing averaged potential as 'eam' (funcfl) format is experimental.")
        # Need to select the *last* element's data (the A-atom)
        avg_idx = -1
        crystal_structures_str = crystal_structures[avg_idx] # Just take the A-atom structure
        atline = f"{int(atomic_numbers[avg_idx])} {float(atomic_masses[avg_idx]):.4f} {float(lattice_parameters[avg_idx]):.4f} {crystal_structures_str}"
        parameterline = f'{int(Nrho)}\t{float(drho):.16e}\t{int(Nr)}\t{float(dr):.16e}\t{float(cutoff):.10e}'
        potheader = f"# EAM potential (Averaged Atom 'X') from : # {source} \n {atline} \n {parameterline}"
        potfile = open(out_file,'wb')
        # Write F (only the last row for A-atom)
        np.savetxt(potfile, F[avg_idx,:].reshape(Nrho), fmt='%.16e', header=potheader, comments='')
        # Write charge Z (calculated from rep[avg_idx, avg_idx, :])
        charge = rep[avg_idx, avg_idx, :] / _hartree_bohr_in_electronvolt_angstrom
        # Add small epsilon to avoid sqrt(0) or sqrt(negative) due to float precision
        charge[charge < 0] = 0
        charge = np.sqrt(charge + 1e-99)
        np.savetxt(potfile, charge.reshape(Nr), fmt='%.16e')
        # Write f (only the last row/component for A-atom - needs care depending on how f was averaged)
        # Assuming f passed is the 2D (N+1, Nr) array here for simplicity if called directly
        f_to_write = f[avg_idx,:] # Get the A-atom density row
        np.savetxt(potfile, f_to_write.reshape(Nr), fmt='%.16e')
        potfile.close()
        # --- End of experimental 'eam' block ---


    elif kind == "eam/alloy":
        # This block should work correctly with N+1 elements
        potheader = f"# Averaged EAM alloy potential from :\n# {source} \n# \n"
        potfile = open(out_file,'wb')
        np.savetxt(
            potfile, elements, fmt="%s", newline=' ',
            header=potheader+str(num_elements), # num_elements is N+1
            footer=f'\n{Nrho}\t{drho:e}\t{Nr}\t{dr:e}\t{cutoff:e}\n',
            comments=''
        )
        # Write F and f tables (f should be 2D (N+1, Nr) here)
        if f.ndim != 2 or f.shape[0] != num_elements:
             raise ValueError(f"For kind 'eam/alloy', 'f' must be 2D with shape ({num_elements}, {Nr}), but got {f.shape}")
        for i in range(num_elements): # Iterates 0 to N
            np.savetxt(
                potfile, np.append(F[i,:], f[i,:]), fmt="%.16e",
                header=f'{int(atomic_numbers[i])}\t{atomic_masses[i]:.4f}\t{lattice_parameters[i]:.4f}\t{crystal_structures[i]}',
                comments=''
            )
        # Write pair interactions tables (rep is 3D (N+1, N+1, Nr))
        # This list comprehension correctly handles the lower triangle for N+1 elements
        [[np.savetxt(potfile,rep[i,j,:],fmt="%.16e") for j in range(num_elements) if j <= i] for i in range(num_elements)]
        potfile.close()

    elif kind == "eam/fs":
        # This block needs to correctly handle f being 3D (N+1, N+1, Nr)
        potheader = f"# Averaged EAM fs potential from :\n# {source} \n# \n"
        potfile = open(out_file,'wb')
        np.savetxt(
            potfile, elements, fmt="%s", newline=' ',
            header=potheader+str(num_elements), # num_elements is N+1
            footer=f'\n{Nrho}\t{drho:e}\t{Nr}\t{dr:e}\t{cutoff:e}\n',
            comments=''
        )
        # Write F and f tables (f should be 3D (N+1, N+1, Nr) here)
        if f.ndim != 3 or f.shape[0] != num_elements or f.shape[1] != num_elements:
             raise ValueError(f"For kind 'eam/fs', 'f' must be 3D with shape ({num_elements}, {num_elements}, {Nr}), but got {f.shape}")
        for i in range(num_elements): # Iterates 0 to N
            # Flatten the f[i, :, :] slice which has shape (N+1, Nr)
            f_slice_flat = f[i,:,:].flatten()
            data_to_write = np.append(F[i,:], f_slice_flat)
            # Check expected length: Nrho + (N+1)*Nr
            expected_len = Nrho + num_elements * Nr
            if len(data_to_write) != expected_len:
                # This check might be slightly off if the interpretation of the fs format writing is complex
                 print(f"Warning: Length mismatch when writing eam/fs element {i}. "
                       f"Expected {expected_len}, got {len(data_to_write)}. "
                       f"Check format specification if issues arise.")

            np.savetxt(
                potfile, data_to_write, fmt="%.16e",
                header=f'{int(atomic_numbers[i])}\t{atomic_masses[i]:.4f}\t{lattice_parameters[i]:.4f}\t{crystal_structures[i]}',
                comments=''
            )
        # Write pair interactions tables (rep is 3D (N+1, N+1, Nr))
        # This list comprehension correctly handles the lower triangle for N+1 elements
        [[np.savetxt(potfile, rep[i,j,:], fmt="%.16e") for j in range(num_elements) if j <= i] for i in range(num_elements)]
        potfile.close()
    else:
        raise ValueError(f"EAM kind '{kind}' is not supported for writing")


def find_potential(elements, potentials_dir):
    """
    Find appropriate EAM potential file based on alloy elements.
    
    Parameters
    ----------
    elements : list of str
        List of element symbols (e.g., ['Fe', 'Ni', 'Cr']).
    potentials_dir : str or Path
        Directory containing potential files.
    
    Returns
    -------
    potential_file : Path
        Path to the selected potential file.
    kind : str
        Kind of potential, e.g. 'eam/alloy', 'eam/fs', or 'lammps/eam'.
    
    Raises
    ------
    FileNotFoundError
        If no matching potential file is found.
    """
    potentials_path = Path(potentials_dir)
    # Collect files matching all supported patterns
    potential_files = list(potentials_path.glob("*.eam.alloy")) + \
                      list(potentials_path.glob("*.eam.fs")) + \
                      list(potentials_path.glob("*.lammps.eam")) + \
                      list(potentials_path.glob("*.eam"))
    
    # Normalize element symbols for comparison
    query_elements = [el.lower() for el in elements]
    query_key = "".join(query_elements)  # e.g., 'fenicr'
    
    exact_matches = []
    substring_matches = []
    
    for file in potential_files:
        # Extract the base part of the filename before the first dot.
        # It is assumed that the filename follows the format FeNiCr.eam.alloy, where "FeNiCr" is the alloy designation.
        base = file.name.split('.')[0].lower()  # e.g., 'fenicr' or 'fenicrmn'
        
        if base == query_key:
            exact_matches.append(file)
        elif all(el in base for el in query_elements):
            substring_matches.append(file)
    
    # Prioritize exact matches if any are found.
    if exact_matches:
        matches = exact_matches
    elif substring_matches:
        matches = substring_matches
    else:
        raise FileNotFoundError(f"No suitable potential found for elements: {', '.join(elements)}")
    
    # If multiple variants are found, select the file with the shortest name,
    # which usually corresponds to the exact or least "overloaded" name.
    matches.sort(key=lambda x: len(x.name))
    eam_potential_path = matches[0]
    
    # Determine potential type based on the file extension
    if eam_potential_path.name.endswith(".eam.alloy"):
        kind = "eam/alloy"
    elif eam_potential_path.name.endswith(".eam.fs"):
        kind = "eam/fs"
    elif eam_potential_path.name.endswith(".lammps.eam"):
        kind = "lammps/eam"
    elif eam_potential_path.name.endswith(".eam"):
        kind = "eam"
    else:
        kind = "unknown"
    
    return eam_potential_path, kind



def average_potential(
    concentrations,
    parameters,
    F,
    f,
    rep,
    kind="eam/alloy", # Keep track of the original kind
    avg_atom="X",     # Changed default symbol to X
    atomic_number=999,
    crystal_structure="unknown",
    lattice_constant=1.0,
):
    r"""Generate Average-atom potential, preserving dimensionality based on input kind."""

    # --- Input Validation --- (Keep the checks from the previous version)
    if kind == "eam" or kind == "lammps/eam":
        raise NotImplementedError(
             f"Cannot create an average-atom potential from a single-element "
             f"potential of kind '{kind}'. Averaging requires a multi-element "
             f"potential (e.g., 'eam/alloy' or 'eam/fs')."
         )
    if not np.isclose(np.sum(concentrations), 1):
        raise ValueError(f"Concentrations must sum to 1. Got: {concentrations} (sum={np.sum(concentrations)})")
    num_original_elements = len(parameters.symbols)
    if len(concentrations) != num_original_elements:
        raise ValueError(f"Number of concentrations ({len(concentrations)}) must match "
                         f"number of elements in the potential ({num_original_elements}).")

    # --- Parameter Setup --- (Same as before)
    symbols = [s for s in parameters.symbols] + [avg_atom]
    atomic_numbers = np.r_[parameters.atomic_numbers, atomic_number]
    atomic_masses = np.r_[
        parameters.atomic_masses, np.average(parameters.atomic_masses, weights=concentrations)
    ]
    # Use calculated average lattice constant unless overridden
    avg_lattice_constant = np.average(parameters.lattice_constants, weights=np.array(concentrations))
    lattice_constants = np.r_[parameters.lattice_constants, lattice_constant if lattice_constant != 1.0 else avg_lattice_constant]
    crystal_structures = np.r_[
        parameters.crystal_structures, np.array(crystal_structure)
    ]
    new_parameters = EAMParameters(
        symbols, atomic_numbers, atomic_masses, lattice_constants, crystal_structures,
        parameters.number_of_density_grid_points, parameters.number_of_distance_grid_points,
        parameters.density_grid_spacing, parameters.distance_grid_spacing, parameters.cutoff,
    )

    # --- Average Embedding Energy F --- (Same as before)
    new_F = np.zeros((num_original_elements + 1, F.shape[1]), dtype=F.dtype)
    new_F[:-1, :] = F
    new_F[-1, :] = np.average(F, axis=0, weights=concentrations)

    # --- Average Electron Density f (Return 2D or 3D based on input kind) ---
    Nr = parameters.number_of_distance_grid_points
    if f.ndim == 2 and kind == "eam/alloy":
        # Return 2D: (N+1, Nr)
        new_f = np.zeros((num_original_elements + 1, Nr), dtype=f.dtype)
        new_f[:-1, :] = f
        new_f[-1, :] = np.average(f, axis=0, weights=concentrations)

    elif f.ndim == 3 and kind == "eam/fs":
        # Return 3D: (N+1, N+1, Nr)
        new_f = np.zeros((num_original_elements + 1, num_original_elements + 1, Nr), dtype=f.dtype)
        new_f[:-1, :-1, :] = f # Copy original N x N block

        # Calculate the single average density function g_A(r)
        avg_f_contrib_X = np.average(f, axis=1, weights=concentrations) # Shape (N, Nr)
        avg_f_A = np.average(avg_f_contrib_X, axis=0, weights=concentrations) # Shape (Nr)

        # Populate the last row and column (index N) with the average density avg_f_A
        # This assumes density at Y from A (f_AY) = density at A from X (f_XA) = density at A from A (f_AA) = avg_f_A
        new_f[num_original_elements, :, :] = avg_f_A  # Row N (A-atom contribution)
        new_f[:, num_original_elements, :] = avg_f_A  # Column N (density experienced by A-atom)
        # Corner element new_f[N, N, :] is set correctly by the above lines

    else:
        raise ValueError(f"Unexpected combination of electron density array dimensions ({f.ndim}) "
                         f"and potential kind ('{kind}')")

    # --- Average the Pair Potential rep --- (Same as before, returns 3D (N+1, N+1, Nr))
    new_rep = np.zeros((num_original_elements + 1, num_original_elements + 1, Nr), dtype=rep.dtype)
    new_rep[:-1, :-1, :] = rep
    new_rep[-1, :-1, :] = np.average(rep, axis=0, weights=concentrations)
    new_rep[:-1, -1, :] = new_rep[-1, :-1, :]
    new_rep[-1, -1, :] = np.average(new_rep[-1, :-1, :], axis=0, weights=concentrations)

    # Return the new parameters and arrays (new_f will be 2D or 3D)
    return new_parameters, new_F, new_f, new_rep

def generate_averaged_potential(concentrations, eam_potential_path, kind, components, alloy_name):
    """Generate averaged potential for given composition, preserving input format."""
    avg_potentials_dir = Path('./data/input/potentials/average_potentials')
    avg_potentials_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading potential file: {eam_potential_path} (kind: {kind})")
    source, parameters, F, f, rep = read_eam(str(eam_potential_path), kind)

    print(f"Generating average potential for composition: {dict(zip(components, concentrations))}")
    concentrations_arr = np.array(concentrations)

    # Determine lattice constant for the average atom (use weighted average if not specified)
    # Note: average_potential now handles this default internally if lattice_constant=1.0
    
    new_parameters, new_F, new_f, new_rep = average_potential(
        concentrations_arr,
        parameters,
        F,
        f,
        rep,
        kind=kind, # Pass the original kind here
        avg_atom="X",
        atomic_number=999,
        crystal_structure="unknown",
        # lattice_constant=1.0 # Use default average calculated inside average_potential
    )

    composition_str = ", ".join(f"{c*100:.1f}% {e}" for e, c in zip(components, concentrations))
    source_info = f"{source}\nAveraged from {kind} potential for composition: {composition_str}"

    # --- Use the original input 'kind' for the output ---
    output_kind = kind
    output_file_suffix = output_kind.replace('/', '.') # e.g., eam-alloy or eam-fs
    output_file = avg_potentials_dir / f"{alloy_name}.averaged.{output_file_suffix}"

    print(f"Writing averaged potential (as {output_kind}) to: {output_file}")
    write_eam(
        source_info,
        new_parameters, # Includes N+1 elements
        new_F,          # Shape (N+1, Nrho)
        new_f,          # Shape (N+1, Nr) or (N+1, N+1, Nr) depending on output_kind
        new_rep,        # Shape (N+1, N+1, Nr)
        str(output_file),
        kind=output_kind # Use the original kind for writing
    )

    print(f"Successfully generated averaged potential: {output_file}")
    return output_file
