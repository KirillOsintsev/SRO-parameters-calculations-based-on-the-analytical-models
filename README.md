# Short-Range Order (SRO) Parameters calculator based on the analytical models

## Overview

This computational framework performs calculation of the Short-Range Order (SRO) parameters, Effective Pair Interactions (EPIs), misfit volumes, and normalized SRO parameters in the metallic systems composed of 2-6 elements. The framework integrates molecular dynamics simulations (via LAMMPS) with analytical models based on the theories of Y. Rao, W. Curtin (DOI: 10.1016/j.actamat.2022.117621) and S. Nag, W. Curtin (DOI: 10.1016/j.actamat.2023.119472) to predict ordering behavior in complex alloys across a range of temperatures.

## Usage
1. Ensure input files exist:
   - `data/input/database.json`
   - `data/input/{name-of-your-excel-file}.xls` # Note: only .xls file format is currently supported
   - EAM potential in `data/input/potentials/eam/`
2. Run the interactive CLI:

Open terminal and run
```
python main.py
```
The command line interface (CLI) will be run. 

### Interacting with the CLI

Once the script starts, you will be prompted to provide the following information step-by-step:

**Operation Mode**: Choose whether to perform a fresh simulation or load cached data.

1. **Calculate new**: Runs the full SRO parameter simulation.

2. **Load from existing files**: Processes data from previously generated results in the ./data/output/ directory.

3. **Alloy System**: Enter the name of the alloy system (e.g., `FeNiCr`).

4. **Crystal Phase**: Select `fcc` or `bcc`.

5. **Nearest Neighbors**: Specify number of shells to calculate (1-4). However, currently is implemented only the 1st NN.

6. **Unit Cells**: Number of unit cells for LAMMPS simulations (default: 5). The more the unit cells the higher the computational cost.

7. **Lattice Constant**: (Mode 1 only) Enter lattice constant value in (Å) which will be used for the molecular dynamics energy minimization for all compositions.

7. **Composition Input**: (Mode 1 only) Provide path to Excel file with compositions. Copy and paste the link to your .xls file in the command line. An example of the .xls file is provided in the ./data/input/FeNiCr.xls. The script supports consequent run for the multiple alloy systems written in the first column of the .xls file.

**Configuration**: Temperature Range

The simulation iterates through a predefined list of temperatures to calculate short-range order parameters which can be manufally changed.

**Default Range**: 300 K to 2000 K (including specific points like 773 K and 1273 K).

**How to Modify**: The temperature points are hardcoded in the main() function within main.py. To change the range or add specific temperature steps, locate the temperatures list and update it:

```
# main.py

def main():
    # Modify this list to change the temperature range
    temperatures = [300, 400, 500, ..., 2000]
```

## Project Structure

```
Project/
├── main.py                   # Main execution script (interactive CLI)
├── calculators/                        # Core calculation modules
│   ├── __init__.py
│   ├── EPI_calculator.py               # Effective Pair Interaction (EPI) calculations
│   ├── calculate_misfit_volume.py      # Misfit volume calculations (from database.json)
│   ├── calculate_sro_parameters.py     # Warren–Cowley SRO parameter calculations
│   └── normalized_sro_parameters.py    # Normalized SRO parameters β_i/β_0
├── utils/                              # Utility and helper modules
│   ├── __init__.py
│   ├── composition.py                  # Composition parsing and naming
│   ├── create_calculation_structure.py # Directory structure creation + CLI helpers
│   ├── generate_calculation_files.py   # LAMMPS/solver script templates (generated per pair/shell)
│   ├── create_average_potential.py     # EAM potential averaging utilities
│   └── normalized_sro_screening.py     # Normalized SRO screening + report generation
└── data/
    ├── input/
    │   ├── database.json          # Element property database
    │   ├── potentials/
    │   │   ├── eam/              # EAM potential files
    │   │   └── average_potentials/ # Generated averaged potentials
    │   └── *.xls                 # Composition input files
    └── output/
        └── results/
            └── {alloy_system}/   # Results organized by alloy system
                └── {alloy_name}/ # Individual alloy results
                    ├── {alloy_name}.json              # Configuration file
                    ├── {alloy_name}_epi_results.txt   # EPI results
                    ├── {alloy_name}_misfit_volumes_results.txt # Misfit volumes
                    |── {alloy_name}_pair_probabilities.json  # Pair probabilities
                    ├── {alloy_name}_sro_results.json  # SRO parameters
                    ├── {alloy_name}_fugacity_results.txt # Fugacity values
                    ├── {alloy_name}_normalized_sro_parameters.txt # Normalized SRO
                    ├── {pair_name}/                   # Pair-specific directories
                    │   └── {nn}/                     # Nearest neighbor shells
                    │       ├── calculate_{nn}_S-S_interaction_by_pressure_relaxation.py
                    │       ├── run_calculation_{nn}.sh
                    │       └── log*.lammps
                    └── all_alloys_normalized_sro_report.txt # Summary report
```

## Methodology

### 1. Composition Processing

The system accepts alloy compositions from the Excel file:
- **Excel file input**: Batch processing from `.xls` files containing composition data

Compositions are parsed generating:
- Alloy name (e.g., `Fe56Ni23Cr21`)
- Component list and concentrations
- Number of components and unique pairs

### 2. Potential Generation

For each composition, the system:
1. **Locates EAM potential**: Searches for appropriate Embedded Atom Method (EAM) potential files matching the alloy components
2. **Generates averaged potential**: Creates a concentration-weighted average potential using the `generate_averaged_potential` function
3. **Saves averaged potential**: Stores the averaged potential in `data/input/potentials/average_potentials/`

### 3. Calculation Structure Setup

The system creates a hierarchical directory structure:
```
{alloy_name}/
├── {pair_name}/          # For each unique element pair (e.g., CrFe, CrNi, FeNi)
    └── {nn}/             # For each nearest neighbor shell (1st, 2nd)
        ├── calculate_{nn}_S-S_interaction_by_pressure_relaxation.py
        └── run_calculation_{nn}.sh
```

### 4. LAMMPS Calculations

For each pair and nearest neighbor shell, the system:
1. **Generates Python scripts**: Creates LAMMPS calculation scripts for solute-solute interactions
2. **Generates shell scripts**: Creates bash scripts to execute calculations for all pair combinations
3. **Executes calculations**: Runs LAMMPS simulations via subprocess calls
4. **Extracts results**: Parses log files to extract interaction energies


### 5. Effective Pair Interactions (EPIs)

The `EPI_calculator` module:
- Parses LAMMPS log files from each pair/shell combination
- Calculates effective pair interaction energies: `V_eff = E_12 - (E_11 + E_22)/2`
- Aggregates results across all pairs and shells
- Saves results to `{alloy_name}_epi_results.txt`

### 6. Misfit Volume Calculation

The `misfit_volume_calculator` module:
- Loads atomic volumes from `database.json`
- Calculates average atomic volume: `V_avg = Σ(c_i × V_i)`
- Computes misfit volumes: `ΔV_i = V_i - V_avg`
- Saves results to `{alloy_name}_misfit_volumes_results.txt`

**Note**: The use of a composition-weighted average volume is consistent with a (Vegard-like) linear mixture rule
commonly used in solute-strengthening models, e.g. Varvenne et al. (2016), DOI: `10.1016/j.actamat.2016.07.040`.

### 7. SRO Parameter Calculation

The `calculate_sro_parameters` module:
- **Fugacity calculation**: Computes temperature-dependent fugacity values from EPIs
- **SRO parameter calculation**: Determines Warren-Cowley SRO parameters (α_nm) for each pair and shell based on the analytical models derived in the Y.Rao and W.Curtin's article (DOI: 10.1016/j.actamat.2022.117621). Currently supports prediction of the SRO parameters only for the 1st coordination shell for the binary, ternary, quaternary and quinary systems using pair models. 
  
- **Temperature dependence**: Calculates SRO parameters across the specified temperature range
- Saves results to `{alloy_name}_sro_results.json` and `{alloy_name}_fugacity_results.txt`

### 8. Normalized SRO Parameters

The `normalized_sro_parameters` module is based on the solute strengthening theory for alloys with SRO developed by S.Nag and W. Curtin DOI: '10.1016/j.actamat.2023.119472':
- **Beta calculation**: Computes normalized SRO parameters (β_i) incorporating misfit volumes:
  ```
  β_i = Σ_n,m c_n × c_m × (ΔV_n - ΔV_m)² × α_nm^(i)
  ```
- **Beta_0 normalization**: Calculates reference value:
  ```
  β_0 = Σ_i c_i × (ΔV_i)²
  ```
- **Phase-specific averaging**:
  - **FCC**: `β_1/β_0` and `(β_2 + 4β_3 + 2β_4)/β_0`
  - **BCC**: `(1.5β_1 + β_2)/β_0` and `(β_3 + 1.8β_4)/β_0`
- Saves results to `{alloy_name}_normalized_sro_parameters.txt`

### 9. Report Generation

The system generates comprehensive reports:
- **Individual alloy reports**: Detailed results for each composition
- **Summary report**: Aggregated results across all alloys in the system (`all_alloys_normalized_sro_report.txt`)
- Reports include normalized SRO parameters, condition satisfaction flags, and temperature-dependent behavior


### Output Files

#### Configuration File (`{alloy_name}.json`)
Contains all calculation parameters:
- Alloy name and composition
- Component concentrations
- Nearest neighbor shells
- Lattice constant and phase
- Potential style

#### EPI Results (`{alloy_name}_epi_results.txt`)
Effective pair interaction energies for each shell and pair combination.

#### Misfit Volumes (`{alloy_name}_misfit_volumes_results.txt`)
Misfit volumes for each element in the alloy.

#### Pair Probabilities (`{alloy_name}_pair_probabilities.json`)
Temperature-dependent pair probabilities for all pairs and shells used in SRO calculations.

#### SRO Parameters (`{alloy_name}_sro_results.json`)
Temperature-dependent Warren-Cowley SRO parameters (α_nm) for all pairs and shells.

#### Fugacity Results (`{alloy_name}_fugacity_results.txt`)
Temperature-dependent fugacity values used in SRO calculations.

#### Normalized SRO Parameters (`{alloy_name}_normalized_sro_parameters.txt`)
Normalized SRO parameters (β_i/β_0) and phase-specific averages.

#### Summary Report (`all_alloys_normalized_sro_report.txt`)
Aggregated results for all alloys in the system.

## Dependencies

### Python Packages
- `numpy`: Numerical computations
- `pandas`: Excel file reading
- `pathlib`: File path handling
- `json`: Configuration and data storage
- `subprocess`: LAMMPS execution
- `matplotlib`: Visualization (optional)
- `scipy`: Optimization routines (SRO solving; also used in generated scripts)
- `mpi4py`: MPI bindings for Python (used by generated LAMMPS scripts)
- `lammps`: Python interface to LAMMPS (needed by generated scripts)

Install the Python dependencies in the same environment used to run this project:

```bash
pip install -r requirements.txt
```

### Python Version
- Recommended: **Python 3.8+** (SciPy/matplotlib compatibility).

### External Software
- **LAMMPS**: Molecular dynamics simulator (https://www.lammps.org/)
- **mpi4py**: Parallel execution support (optional, used by generated scripts)

### Data Requirements
- **EAM Potentials**: Embedded Atom Method potential files in LAMMPS format
- **Database File**: `database.json` containing element properties (atomic volumes, lattice constants, etc.)

## Key Functions

### Main Functions (`main.py`)

- `process_composition()`: Orchestrates the complete calculation pipeline for a single composition
- `main()`: User interface and batch processing coordinator

### Calculator Modules (`calculators/`)

- `calculate_epis()`: Extracts and calculates effective pair interactions
- `misfit_volume_calculator()`: Computes misfit volumes from atomic volumes
- `calculate_sro_parameters()`: Determines Warren-Cowley SRO parameters
- `calculate_normalized_sro_parameters()`: Calculates normalized SRO parameters

### Utility Modules (`utils/`)

- `get_alloy_name_and_concentrations()`: Parses composition and generates naming
- `create_folder_structure()`: Sets up calculation directory hierarchy
- `find_potential()`: Locates appropriate EAM potential files
- `generate_averaged_potential()`: Creates concentration-weighted averaged potentials
- `create_shell_script()` / `create_python_script()`: Generates LAMMPS calculation scripts
- `save_normalized_sro_parameters()` / `load_normalized_sro_parameters()`: I/O for normalized SRO data
- `generate_alloy_report()`: Creates summary reports

## Scientific Background

### Short-Range Order (SRO)

Warren-Cowley SRO parameters (α_nm) quantify the deviation from random mixing:
- α_nm = 0: Random distribution
- α_nm > 0: Clustering tendency
- α_nm < 0: Ordering tendency

### Normalized SRO Parameters

Normalized parameters (β_i/β_0) incorporate chemical and size effects:
- **β_i**: Shell-specific normalized SRO incorporating misfit volumes
- **β_0**: Reference normalization factor
- Phase-specific averaging accounts for crystal structure effects

**Literature basis**:
- Rao, Y., Curtin, W. A. (2022), *Analytical models of short-range order in FCC and BCC alloys*, Acta Materialia 226, 117621. DOI: `10.1016/j.actamat.2022.117621`
- Nag, S. and Curtin, W. A. (2024), *Solute-strengthening in metal alloys with short-range order*, Acta Materialia 263, 119472. DOI: `10.1016/j.actamat.2023.119472`

### Effective Pair Interactions

EPIs (V_nm) represent the energy difference between mixed and separated pairs:
- Directly related to ordering/clustering tendencies
- Temperature-dependent through fugacity relationships
- Used to predict SRO behavior via statistical thermodynamics

## Notes

- The system supports binary, ternary, quaternary, and quinary alloys
- Results are cached; existing calculation files are reused when available
- The system automatically handles file existence checks to avoid redundant calculations

## Authors & Contributions (code)

- **Osintsev Kirill (Wenzhou university)**
  - Implemented core pipeline modules.
- **Xin Liu (EPFL)** — `xin.liu@epfl.ch`
  - Implemented run-script template used by `utils/generate_calculation_files.py`.

## Third-party Code & Licenses

- `utils/create_average_potential.py` includes code adapted from `matscipy` (`average_atom.py`):
  - Lars Pastewka (U. Freiburg), 2021
  - Wolfram G. Nöhring (U. Freiburg), 2020
  - Project: `matscipy` — `https://github.com/libAtoms/matscipy`
  - License: GNU General Public License (GPL), v2 or later (see `http://www.gnu.org/licenses/`)

## References

This implementation follows established methods for:
- Warren-Cowley SRO parameter calculation
- EAM potential averaging techniques
- Statistical thermodynamics of ordering in alloys
- Misfit volume calculations for multi-component systems

Key papers referenced by the implementation:
- Rao, Y., Curtin, W. A. (2022), *Analytical models of short-range order in FCC and BCC alloys*, Acta Materialia 226, 117621. DOI: `10.1016/j.actamat.2022.117621`
- Nag, S. and Curtin, W. A. (2024), *Solute-strengthening in metal alloys with short-range order*, Acta Materialia 263, 119472. DOI: `10.1016/j.actamat.2023.119472`
- Varvenne, C., Luque, A., Curtin, W. A. (2016), *Theory of strengthening in fcc high entropy alloys*, Acta Materialia 118, 164–176. DOI: `10.1016/j.actamat.2016.07.040`
