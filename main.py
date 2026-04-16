import utils
import calculators
from pathlib import Path
import subprocess
import json

def process_composition(composition, database_path, lattice_constant, temperatures, phase, nn_list, ncells, alloy_system, recalculate=True):
    """
    Process a composition, either calculating or loading normalized SRO parameters.
    
    Args:
        composition (dict): {element: fraction}
        database_path (Path): Path to database file
        temperatures (list): Temperature range
        phase (str): 'fcc' or 'bcc'
        selected_model (str): Strength model
        nn_list (list): Nearest neighbor shells
        ncells (int): Number of unit cells
        alloy_system (str): Alloy system folder
        recalculate (bool): True to calculate, False to load existing data
    
    Returns:
        dict: Alloy data including normalized SRO parameters
    """
    alloy_name, concentrations = utils.get_alloy_name_and_concentrations(composition)
    components = list(composition.keys())
    number_of_components = len(components)
    
    alloy_path = Path(f'./data/output/results/{alloy_system}/{alloy_name}')
    alloy_path.mkdir(parents=True, exist_ok=True)
    
    if recalculate:
        # Lattice constant
       # lattice_constant = utils.get_lattice_constant(components)
        
        # Create pairs and folder structure
        pairs = utils.create_folder_structure(alloy_path, components, nn_list)
        
        # Generate average potential
        potentials_dir = Path('./data/input/potentials/eam')
        eam_potential_path, kind = utils.find_potential(components, potentials_dir)
        utils.generate_averaged_potential(concentrations, eam_potential_path, kind, components, alloy_name)
        
        # Configuration
        config = {
            "alloy_name": alloy_name, "composition": composition, "concentrations": concentrations,
            "components": components, "number_of_components": number_of_components,
            "pairs": [list(pair) for pair in pairs], "nearest_neighbors": nn_list,
            "lattice_constant": lattice_constant, "ncells": ncells, "phase": phase,
            "potential_style": kind
        }
        with open(alloy_path / f'{alloy_name}.json', 'w') as f:
            json.dump(config, f, indent=4)
        
        # EPI calculation
        epi_file = alloy_path / f'{alloy_name}_epi_results.txt'
        if not epi_file.exists():
            calculation_scripts = []
            for pair in pairs:
                pair_name = ''.join(pair)
                pair_path = alloy_path / pair_name
                for nn in nn_list:
                    script = utils.create_shell_script(pair_path, nn, config, pair)
                    calculation_scripts.append(script)
                    utils.create_python_script(pair_path, nn, config)
            for script in calculation_scripts:
                result_file = script.parent / 'log.lammps'
                if not result_file.exists():
                    subprocess.run(['bash', script.name], cwd=script.parent, check=True, capture_output=True, text=True)
            calculators.calculate_epis(alloy_path, alloy_name, nn_list)
        
        # Misfit volumes
        misfit_volume_file = alloy_path / f'{alloy_name}_misfit_volumes_results.txt'
        if not misfit_volume_file.exists():
            calculators.misfit_volume_calculator(database_path, alloy_path, alloy_name, composition, components, concentrations)
        
        # SRO parameters
        sro_output_file = alloy_path / f'{alloy_name}_sro_results.json'
        fugacity_output_file = alloy_path / f'{alloy_name}_fugacity_results.txt'
        if not sro_output_file.exists():
            sro_parameters = calculators.calculate_sro_parameters(
                epi_file, sro_output_file, concentrations, number_of_components, fugacity_output_file,
                nn_list, components, alloy_name, temperatures, alloy_path
            )
        else:
            with open(sro_output_file, 'r') as f:
                sro_parameters = json.load(f)
        
        # Normalized SRO parameters
        normalized_sro_file = alloy_path / f'{alloy_name}_normalized_sro_parameters.txt'
        if not normalized_sro_file.exists():
            averages = calculators.calculate_normalized_sro_parameters(
                sro_parameters, sro_output_file, misfit_volume_file, phase, composition, alloy_path, alloy_name
            )
            utils.save_normalized_sro_parameters(averages, normalized_sro_file)
    
    # Load normalized SRO parameters
    normalized_sro_file = alloy_path / f'{alloy_name}_normalized_sro_parameters.txt'
    averages = utils.load_normalized_sro_parameters(normalized_sro_file) if normalized_sro_file.exists() else {}
    
    return {
        'alloy_name': alloy_name,
        'composition': composition,
        'averages': averages
    }

def main():
    database_path = Path('./data/input/database.json')
    #temperatures = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
    temperatures = [300, 400, 500, 600, 700, 773, 800, 900, 1000, 1100, 1200, 1273, 1400, 1500, 1600, 1700, 1800, 1900, 2000]
    #temperatures = [350, 450, 550, 650, 750, 850, 950, 1050, 1150, 1250]
    #temperatures = [1000]
    results = []
    
    mode = input("Choose mode: (1) Calculate new, (2) Load from existing files: ").strip()
    recalculate = mode == '1'
    
    alloy_system = input("Enter alloy system: ").strip()
    phase = utils.get_phase()
    nn_list = utils.get_user_input_nearest_neighbors()
    ncells = utils.get_ncells()
    """
    print("Available models:")
    for i, model in enumerate(calculators.single_calc_strength().ssmodels_all, start=1):
        print(f"{i}: {model}")
    model_choice = input("Select a model by number: ")
    try:
        selected_model = calculators.single_calc_strength().ssmodels_all[int(model_choice) - 1]
    except (ValueError, IndexError):
        print("Invalid choice, defaulting to 'FCC_Varvenne-Curtin-2016'.")
        selected_model = calculators.single_calc_strength().ssmodels_all[0]
    """
    if recalculate:
        
        lattice_constant = float(input("Enter the lattice constant to use for all compositions: "))
        compositions = utils.load_compositions_from_excel()
        for composition in compositions:
            result = process_composition(composition, database_path, lattice_constant, temperatures, phase,
                                        nn_list, ncells, alloy_system, recalculate)
            if result:
                results.append(result)
    else:
        base_path = Path(f'./data/output/results/{alloy_system}')
        if not base_path.exists():
            print(f"Folder {base_path} not found.")
            return
        for alloy_dir in base_path.iterdir():
            if alloy_dir.is_dir():
                alloy_name = alloy_dir.name
                config_file = alloy_dir / f'{alloy_name}.json'
                if config_file.exists():
                    with open(config_file, 'r') as f:
                        config = json.load(f)
                    composition = config['composition']
                    lattice_constant = config.get('lattice_constant', 3.5)  # Use default if not in config
                    result = process_composition(composition, database_path, lattice_constant, temperatures, phase,
                                                nn_list, ncells, alloy_system, recalculate)
                    if result:
                        results.append(result)
    
    if results:
        output_path = Path(f'./data/output/results/{alloy_system}')
        utils.generate_alloy_report(results, output_path, phase, temperatures)

if __name__ == "__main__":
    main()
