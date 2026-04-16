import os
import json
from pathlib import Path

# Existing functions assumed: get_alloy_name_and_concentrations, create_folder_structure, find_potential,
# generate_averaged_potential, get_lattice_constant, create_shell_script, create_python_script,
# get_user_input_nearest_neighbors, get_ncells, get_phase, get_composition_from_user, load_compositions_from_excel

def save_normalized_sro_parameters(averages, output_file):
    """Save normalized SRO parameters to a text file."""
    with open(output_file, 'w') as f:
        for temp, metrics in averages.items():
            f.write(f"Temperature {temp} K\n")
            for key, value in metrics.items():
                f.write(f"{key}: {value}\n")
            f.write("\n")

def load_normalized_sro_parameters(file_path):
    """Load normalized SRO parameters from a text file."""
    averages = {}
    if not file_path.exists():
        return averages
    with open(file_path, 'r') as f:
        lines = f.readlines()
    current_temp = None
    for line in lines:
        line = line.strip()
        if line.startswith("Temperature"):
            current_temp = line.split()[1]
            averages[current_temp] = {}
        elif line and current_temp and ':' in line:
            key, value = line.split(': ', 1)
            averages[current_temp][key] = float(value)
    return averages

def check_normalized_sro_condition(averages, phase):
    """Check if beta_1/beta_0 >= 0.025 and average beta/beta_0 > 0 at any temperature."""
    if phase == "fcc":
        for temp, metrics in averages.items():
            if ("beta_1/beta_0" in metrics and "beta_average/beta_0" in metrics and
                metrics["beta_1/beta_0"] >= 0.025 and metrics["beta_average/beta_0"] > 0):
                return True
    elif phase == "bcc":
        for temp, metrics in averages.items():
            if ("beta_1_average/beta_0" in metrics and "beta_2_average/beta_0" in metrics and
                metrics["beta_1_average/beta_0"] >= 0.05 and metrics["beta_2_average/beta_0"] > 0):
                return True
    return False

def generate_alloy_report(results, output_path, phase, temperatures):
    """Generate a single report for all alloys with normalized SRO parameters."""
    report_file = output_path / 'all_alloys_normalized_sro_report.txt'
    with open(report_file, 'w') as f:
        f.write(f"All Alloys Normalized SRO Report\n")
        f.write(f"Alloy System: {output_path.name}\n")
        f.write(f"Phase: {phase}\n")
        f.write(f"Temperature Range: {min(temperatures)} K to {max(temperatures)} K\n")
        if phase == "fcc":
            f.write(f"Condition: beta_1/beta_0 >= 0.025 and beta_average/beta_0 > 0\n\n")
        elif phase == "bcc":
            f.write(f"Condition: beta_1_average/beta_0 >= 0.05 and beta_2_average/beta_0 > 0\n\n")
        
        for result in sorted(results, key=lambda x: x['alloy_name']):
            alloy_name = result['alloy_name']
            composition = result['composition']
            averages = result['averages']
            satisfies_condition = check_normalized_sro_condition(averages, phase)
            
            f.write(f"Alloy: {alloy_name}\n")
            f.write(f"Composition: {composition}\n")
            f.write("Normalized SRO Parameters:\n")
            for temp in temperatures:
                temp_str = str(temp)
                if temp_str in averages:
                    metrics = averages[temp_str]
                    if phase == "fcc":
                        beta_avg_key = "beta_average/beta_0"
                        beta1_key = "beta_1/beta_0"
                        beta1 = metrics.get(beta1_key, 0.0)
                        beta_avg = metrics.get(beta_avg_key, 0.0)
                        condition_met = beta1 >= 0.025 and beta_avg > 0
                    elif phase == "bcc":
                        beta1_key = "beta_1_average/beta_0"
                        beta_avg_key = "beta_2_average/beta_0"
                        beta1 = metrics.get(beta1_key, 0.0)
                        beta_avg = metrics.get(beta_avg_key, 0.0)
                        condition_met = beta1 >= 0.05 and beta_avg > 0
                    
                    beta1_str = f"_{beta1:.17f}_" if condition_met and satisfies_condition else f"{beta1:.17f}"
                    beta_avg_str = f"_{beta_avg:.17f}_" if condition_met and satisfies_condition else f"{beta_avg:.17f}"
                    f.write(f"  Temperature {temp} K:\n")
                    f.write(f"    {beta_avg_key}: {beta_avg_str}\n")
                    f.write(f"    {beta1_key}: {beta1_str}\n")
            f.write("\n")
    print(f"Report generated at: {report_file}")

# Placeholder if not present
def get_lattice_constant_auto():
    """Placeholder for automatic lattice constant calculation."""
    return 3.3  # Replace with actual logic