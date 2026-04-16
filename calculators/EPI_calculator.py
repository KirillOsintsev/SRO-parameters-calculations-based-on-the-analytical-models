import os
from pathlib import Path
from typing import List, Tuple, Optional

def calculate_epis(alloy_path: str, alloy_name: str = None, 
                   nn_list: List[str] = None) -> List[Tuple[str, str, float, float, float, float]]:

    """
    Calculate EPIs from existing calculation results.
    
    Args:
        base_output_path (str): Base path for output directories.
        alloy_name (str, optional): Name of the alloy.
        logs_dir (str, optional): Direct path to the directory containing calculation results.
        nn_list (List[str], optional): List of nearest neighbors.
        
    Returns:
        List[Tuple[str, str, float, float, float, float]]: Results of EPI calculations.
    """

    
    epi_path = Path(alloy_path)
   
    if not epi_path.exists():
        raise FileNotFoundError(f"Calculation directory not found: {epi_path}")

    output_dir = epi_path
    output_dir.mkdir(exist_ok=True)

    all_results = []
    
    # Process each pair directory
    for pair_dir in epi_path.iterdir():
        if pair_dir.is_dir():
            pair_name = pair_dir.name
            print(f"\nProcessing EPIs for pair: {pair_name}")
            
            # Call automate_epi_calculation
            pair_results = automate_epi_calculation(
                pair_name=str(pair_name),
                logs_dir=str(pair_dir),  
                output_dir=str(output_dir / pair_name),
                nn_list=nn_list  
            )
            all_results.extend(pair_results)

    # Save combined results
    epi_file = save_combined_results(alloy_name, all_results, str(output_dir))
    return all_results, epi_file

def automate_epi_calculation(pair_name: str, logs_dir: str, output_dir: str, nn_list: List[str]) -> List[Tuple[str, str, float, float, float, float]]:
    """
    Calculate EPIs for a specific pair directory.
    
    Args:
        logs_dir (str): Directory containing the log files.
        output_dir (str): Directory to save results.
        nn_list (List[str]): List of nearest neighbors.
        
    Returns:
        List[Tuple[str, str, float, float, float, float]]: Results for this pair.
    """
    results = []

    os.makedirs(output_dir, exist_ok=True)

    for neighbor in nn_list:
        neighbor_path = os.path.join(logs_dir, neighbor)
        log11 = os.path.join(neighbor_path, "log11.txt")
        log22 = os.path.join(neighbor_path, "log22.txt")
        log12 = os.path.join(neighbor_path, "log12.txt")

        # Extract energies
        exx = extract_interaction_energy(log11)
        eyy = extract_interaction_energy(log22)
        exy = extract_interaction_energy(log12)

        if all(e is not None for e in [exx, eyy, exy]):
            veff = calculate_veff(exx, eyy, exy)
            results.append((pair_name, neighbor, "Success", exx, eyy, exy, veff))
           
        else:
            results.append((pair_name, neighbor, "Error: Missing data", 0.0, 0.0, 0.0, 0.0))

    save_pair_results(results, output_dir)

    return results

def extract_interaction_energy(log_file: str) -> Optional[float]:
    """
    Extract interaction energy from a log file.
    
    Args:
        log_file (str): Path to the log file.
        
    Returns:
        Optional[float]: Extracted energy value or None if extraction fails.
    """
    try:
        with open(log_file, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if "#########" in line:
                    energy = float(lines[lines.index(line) + 2].strip())
                    return energy
    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading {log_file}: {e}")
        return None

def calculate_veff(exx: float, eyy: float, exy: float) -> float:
    """
    Calculate effective pair interaction energy.
    
    Args:
        exx (float): Energy for X-X interaction.
        eyy (float): Energy for Y-Y interaction.
        exy (float): Energy for X-Y interaction.
        
    Returns:
        float: Effective pair interaction energy.
    """
    return exx + eyy - 2 * exy

def save_pair_results(results: List[Tuple[str]], output_dir: str):
    """
    Save results for a specific pair to a file.
    
    Args:
        results (List[Tuple[str]]): Results to save.
        output_dir (str): Directory to save results.
    """
    output_file = os.path.join(output_dir, "epi_results.txt")
    
    with open(output_file, 'w') as f:
        f.write("Pair\tNearest Neighbor\tStatus\tExx\tEyy\tExy\tVeff\n")
        
        for neighbor_data in results:
            f.write("\t".join(map(str, neighbor_data)) + "\n")
    
    
    print(f"Pair results saved to {output_file}")

def save_combined_results(alloy_name, all_results: List[Tuple], output_dir: str):
    """
    Save combined results for all pairs.
    
    Args: 
        all_results (List[Tuple]): Combined results for all pairs.
        output_dir (str): Directory to save results.
    """
    output_file = os.path.join(output_dir, f"{alloy_name}_epi_results.txt")
    
    
    with open(output_file, 'w') as f:
        f.write("Pair\tNearest Neighbor\tStatus\tExx\tEyy\tExy\tVeff\n")
        
        for result in all_results:
            f.write("\t".join(map(str, result)) + "\n")
    
    print(f"\nCombined results saved to {output_file}")
    return output_file


if __name__ == "__main__":
    base_output_path = './data/outputs'  
    alloy_name = 'example_alloy'  
    logs_dir = './data/logs/example_alloy'  
    
    try:
        results = calculate_epis(base_output_path=base_output_path,
                                  alloy_name=alloy_name,
                                  logs_dir=logs_dir)
        
        print("EPI Calculation Results:", results)
        
    except Exception as e:
        print(f"An error occurred: {e}")
