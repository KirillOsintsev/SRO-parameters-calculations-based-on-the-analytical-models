import numpy as np
from scipy.optimize import fsolve, root
import matplotlib.pyplot as plt
import json
import argparse
import sys
from pathlib import Path

# --- Constants ---
K_B = 8.617333262145e-5  # Boltzmann constant in eV/K


#Note: the ternary square model is not implemented yet. Only the pair model is implemented.
# Ternary Square Model Constants
SQUARE_MODEL_N_STAR = 4
SQUARE_MODEL_Z = {1: 2, 2: 1}
SQUARE_MODEL_G = [1, 1, 1, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 8, 4, 8, 4, 8]

# Configuration data for atom pairs (distance, (frozenset({atom_indices}), config_index)): count
SQUARE_MODEL_N_IJ_K = {
    1: {
        (frozenset({1, 1}), 1): 4, (frozenset({2, 2}), 2): 4, (frozenset({3, 3}), 3): 4,
        (frozenset({1, 1}), 4): 2, (frozenset({1, 2}), 4): 2, (frozenset({1, 1}), 5): 1,
        (frozenset({1, 2}), 5): 2, (frozenset({2, 2}), 5): 1, (frozenset({1, 2}), 6): 4,
        (frozenset({1, 2}), 7): 2, (frozenset({2, 2}), 7): 2, (frozenset({1, 1}), 8): 2,
        (frozenset({1, 3}), 8): 2, (frozenset({1, 1}), 9): 1, (frozenset({1, 3}), 9): 2,
        (frozenset({3, 3}), 9): 1, (frozenset({1, 3}), 10): 4, (frozenset({1, 3}), 11): 2,
        (frozenset({3, 3}), 11): 2, (frozenset({2, 2}), 12): 2, (frozenset({2, 3}), 12): 2,
        (frozenset({2, 2}), 13): 1, (frozenset({2, 3}), 13): 2, (frozenset({3, 3}), 13): 1,
        (frozenset({2, 3}), 14): 4, (frozenset({2, 3}), 15): 2, (frozenset({3, 3}), 15): 2,
        (frozenset({1, 2}), 16): 2, (frozenset({1, 3}), 16): 2, (frozenset({1, 1}), 17): 1,
        (frozenset({1, 2}), 17): 1, (frozenset({1, 3}), 17): 1, (frozenset({2, 3}), 17): 1,
        (frozenset({1, 2}), 18): 2, (frozenset({2, 3}), 18): 2, (frozenset({1, 2}), 19): 1,
        (frozenset({1, 3}), 19): 1, (frozenset({2, 2}), 19): 1, (frozenset({2, 3}), 19): 1,
        (frozenset({1, 3}), 20): 2, (frozenset({2, 3}), 20): 2, (frozenset({1, 2}), 21): 1,
        (frozenset({1, 3}), 21): 1, (frozenset({2, 3}), 21): 1, (frozenset({3, 3}), 21): 1,
    },
    2: {
        (frozenset({1, 1}), 1): 2, (frozenset({2, 2}), 2): 2, (frozenset({3, 3}), 3): 2,
        (frozenset({1, 1}), 4): 1, (frozenset({1, 2}), 4): 1, (frozenset({1, 2}), 5): 2,
        (frozenset({1, 1}), 6): 1, (frozenset({2, 2}), 6): 1, (frozenset({1, 2}), 7): 1,
        (frozenset({2, 2}), 7): 1, (frozenset({1, 1}), 8): 1, (frozenset({1, 3}), 8): 1,
        (frozenset({1, 3}), 9): 2, (frozenset({1, 1}), 10): 1, (frozenset({3, 3}), 10): 1,
        (frozenset({1, 3}), 11): 1, (frozenset({3, 3}), 11): 1, (frozenset({2, 2}), 12): 1,
        (frozenset({2, 3}), 12): 1, (frozenset({2, 3}), 13): 2, (frozenset({2, 2}), 14): 1,
        (frozenset({3, 3}), 14): 1, (frozenset({2, 3}), 15): 1, (frozenset({3, 3}), 15): 1,
        (frozenset({1, 1}), 16): 1, (frozenset({2, 3}), 16): 1, (frozenset({1, 2}), 17): 1,
        (frozenset({1, 3}), 17): 1, (frozenset({1, 3}), 18): 1, (frozenset({2, 2}), 18): 1,
        (frozenset({1, 2}), 19): 1, (frozenset({2, 3}), 19): 1, (frozenset({1, 2}), 20): 1,
        (frozenset({3, 3}), 20): 1, (frozenset({1, 3}), 21): 1, (frozenset({2, 3}), 21): 1,
    }
}

# --- Parsing and Helper Functions ---

def parse_epi_file(epi_file, nn_list):
    """Parses the EPI results file for Veff values corresponding to nearest neighbors."""
    Veff_values = {distance: {} for distance in nn_list}
    pair_names = []

    try:
        with open(epi_file, 'r') as f:
            lines = f.readlines()[1:]  # Skip header
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                
                pair_name = parts[0]
                if pair_name not in pair_names:
                    pair_names.append(pair_name)
                    
                neighbor = parts[1]
                try:
                    epi_value = float(parts[-1])
                    if neighbor in Veff_values:
                        Veff_values[neighbor][pair_name] = epi_value
                except ValueError:
                    print(f"Warning: Invalid EPI value on line: {line}")
                    continue

    except FileNotFoundError:
        print(f"Error: EPI file not found at {epi_file}")
        sys.exit(1)

    return Veff_values, pair_names


def calculate_fugacity(results_dict, temperatures):
    """Calculates fugacities (x = exp(-Veff / kBT)) for provided temperatures."""
    fugacity_results = {}

    for T in temperatures:
        fugacity_results[T] = {}
        for neighbor, pairs in results_dict.items():
            fugacity_results[T][neighbor] = {}
            for pair_name, Veff in pairs.items():
                fugacity_results[T][neighbor][pair_name] = np.exp(-Veff / (K_B * T))
        
        # Verbose output can be enabled if needed, kept concise for now
        # print(f"Calculated fugacities for {T} K")
    
    return fugacity_results


# --- SRO Calculation Models ---

def binary_sro_pair(concentrations, fugacity_results, pair_names):
    """Calculates SRO for binary alloys using the pair model."""
    c1, c2 = concentrations
    alpha_values = {}
    
    for T, pairs in fugacity_results.items():
        alpha_values[T] = {}
        for neighbor, data in pairs.items():
            x12 = data.get(pair_names[0])

            if x12 is None:
                print(f"Warning: Invalid fugacities for {pair_names[0]} at distance {neighbor}, T={T}")
                continue

            alpha_12 = 1 - (np.sqrt(1 + 4 * c1 * c2 * (x12 - 1)) - 1) / (2 * c1 * c2 * (x12 - 1))
            alpha_values[T][neighbor] = {f'alpha_{pair_names[0]}': alpha_12}
            
    return alpha_values


def ternary_sro_pair(fugacity_results, concentrations, components, probabilities_output_file=None):
    """
    Generic solver for Ternary SRO parameters using the pair approximation.
    Solves mass balance and probability equations dynamically based on component inputs.
    """
    if len(concentrations) != 3 or len(components) != 3:
        raise ValueError("ternary_sro_pair requires exactly 3 concentrations and 3 components.")

    alpha_values = {}
    all_probabilities = {}
    element_to_index = {comp: i for i, comp in enumerate(components)}
    
    def equations(p_vals, x_mat_current, c_arr, comp_indices):
        f = np.zeros(6)
        c0, c1_arg, c2_arg = c_arr[comp_indices[0]], c_arr[comp_indices[1]], c_arr[comp_indices[2]]
        
        # Retrieve fugacities (x01, x02, x12)
        x01 = x_mat_current[comp_indices[0], comp_indices[1]]
        x02 = x_mat_current[comp_indices[0], comp_indices[2]]
        x12 = x_mat_current[comp_indices[1], comp_indices[2]]

        # Probability relation equations
        f[0] = (p_vals[0] * p_vals[3]) / (p_vals[1]**2) - 1/(4 * x01) 
        f[1] = (p_vals[0] * p_vals[5]) / (p_vals[2]**2) - 1/(4 * x02) 
        f[2] = (p_vals[3] * p_vals[5]) / (p_vals[4]**2) - 1/(4 * x12) 

        # Mass Balance Equations
        f[3] = 2 * p_vals[0] + p_vals[1] + p_vals[2] - (2 * c0)
        f[4] = p_vals[1] + 2 * p_vals[3] + p_vals[4] - (2 * c1_arg)
        f[5] = p_vals[2] + p_vals[4] + 2 * p_vals[5] - (2 * c2_arg)
        return f

    for T, neighbors_data in fugacity_results.items():
        alpha_values[T] = {}
        all_probabilities[T] = {}
        
        for neighbor, data in neighbors_data.items():
            # Build fugacity matrix
            x_mat = np.zeros((3, 3))
            missing_fugacity = False
            
            for i in range(3):
                for j in range(i + 1, 3):
                    comp1, comp2 = components[i], components[j]
                    pair_fwd = comp1 + comp2
                    pair_rev = comp2 + comp1
                    
                    val = data.get(pair_fwd, data.get(pair_rev))
                    if val is None:
                        missing_fugacity = True
                        break
                    x_mat[element_to_index[comp1], element_to_index[comp2]] = val
                    x_mat[element_to_index[comp2], element_to_index[comp1]] = val
                if missing_fugacity: break
            
            if missing_fugacity:
                print(f"Skipping {neighbor} at {T} K due to missing fugacity data.")
                continue

            # Solve system
            initial_guess = np.ones(6) / 6
            args = (x_mat, concentrations, tuple(element_to_index[c] for c in components))
            solution, info, ier, msg = fsolve(equations, initial_guess, args=args, full_output=True)

            if ier != 1:
                print(f"Solver failed for {neighbor} at {T} K: {msg}")
                continue

            # Unpack probabilities: p_AA, p_AB, p_AC, p_BB, p_BC, p_CC
            p00, p01, p02, p11, p12, p22 = solution

            all_probabilities[T][neighbor] = {
                f"p_{components[0]}{components[0]}": p00,
                f"p_{components[0]}{components[1]}": p01,
                f"p_{components[0]}{components[2]}": p02,
                f"p_{components[1]}{components[1]}": p11,
                f"p_{components[1]}{components[2]}": p12,
                f"p_{components[2]}{components[2]}": p22
            }
            
            # Calculate Alphas
            a01 = 1 - p01 / (2 * concentrations[0] * concentrations[1])
            a02 = 1 - p02 / (2 * concentrations[0] * concentrations[2])
            a12 = 1 - p12 / (2 * concentrations[1] * concentrations[2])
            
            alpha_values[T][neighbor] = {
                f'alpha_{components[0]}{components[1]}': a01,
                f'alpha_{components[0]}{components[2]}': a02,
                f'alpha_{components[1]}{components[2]}': a12,
            }

            if probabilities_output_file:
                try:
                    with open(probabilities_output_file, 'w') as f:
                        json.dump(all_probabilities, f, indent=4)
                        #print(f"Probabilities saved to {probabilities_output_file}")
                except IOError as e:
                    print(f"Error writing probabilities file: {e}")
            
            # Optional: Print result summary
            #print(f"T={T}K, {neighbor}: alpha_{components[0]}{components[1]}={a01:.4f}, "
                  #f"alpha_{components[0]}{components[2]}={a02:.4f}, alpha_{components[1]}{components[2]}={a12:.4f}")

    return alpha_values


def ternary_sro_square(fugacity_results, concentrations, components, temperatures=None):
    """
    Calculates SRO parameters for a ternary system using the Square Model approximation.
    Requires 1st and 2nd neighbor data.
    """
    element_to_index = {element: idx for idx, element in enumerate(components)}
    
    # Coefficients for balance equations
    coeffs_c1 = [4, 0, 0, 3, 2, 2, 1, 3, 2, 2, 1, 0, 0, 0, 0, 2, 2, 1, 1, 1, 1]
    coeffs_c2 = [0, 4, 0, 1, 2, 2, 3, 0, 0, 0, 0, 3, 2, 2, 1, 1, 1, 2, 2, 1, 1]
    coeffs_c3 = [0, 0, 4, 0, 0, 0, 0, 1, 2, 2, 3, 1, 2, 2, 3, 1, 1, 1, 1, 2, 2]
    
    c1, c2, c3 = concentrations
    alpha_values = {}
    
    for T, neighbors in fugacity_results.items():
        if temperatures and T not in temperatures:
            continue
            
        current_temp_alpha = {} 
        neighbor_map = {"1st": 1, "2nd": 2}
        
        if "1st" not in neighbors or "2nd" not in neighbors:
            print(f"Skipping T={T}K: Missing required neighbor data.")
            continue
        
        # Parse fugacities
        x_values = {1: {}, 2: {}}
        missing_data = False
        
        for neighbor_name, neighbor_data in neighbors.items():
            if neighbor_name not in neighbor_map: continue
            r = neighbor_map[neighbor_name]
            
            for pair, value in neighbor_data.items():
                el1, el2 = None, None
                for c in components:
                    if pair.startswith(c):
                        rem = pair[len(c):]
                        if rem in components:
                            el1, el2 = c, rem
                            break
                if not el1:
                    for c in components:
                        if pair.endswith(c):
                            rem = pair[:-len(c)]
                            if rem in components:
                                el1, el2 = rem, c
                                break
                
                if el1 and el2:
                    idx1 = element_to_index[el1] + 1
                    idx2 = element_to_index[el2] + 1
                    x_values[r][frozenset({idx1, idx2})] = value

        # Validate pairs
        req_pairs = [frozenset({1, 2}), frozenset({1, 3}), frozenset({2, 3})]
        for r in [1, 2]:
            for p in req_pairs:
                if p not in x_values[r]:
                    missing_data = True
        
        if missing_data:
            print(f"Skipping T={T}K: Missing pair fugacities.")
            continue
            
        # Extract fugacities
        x12_1 = x_values[1][frozenset({1, 2})]
        x13_1 = x_values[1][frozenset({1, 3})]
        x23_1 = x_values[1][frozenset({2, 3})]
        x12_2 = x_values[2][frozenset({1, 2})]
        x13_2 = x_values[2][frozenset({1, 3})]
        x23_2 = x_values[2][frozenset({2, 3})]

        # Calculate Energy Term X
        X = 17 * np.log(x12_1 * x13_1 * x23_1) + 7 * np.log(x12_2 * x13_2 * x23_2)
        
        # --- Numerical Solution ---
        def equations(vars):
            p = vars[:21]
            lambda1, lambda2, lambda3 = vars[21:]
            eqs = []
            
            # 21 configuration equations with protection
            for k in range(21):
                p_k = max(p[k], 1e-15)  # Protect against log(0)
                term = (2 * (1 + np.log(p_k)) - 2 * np.log(SQUARE_MODEL_G[k]) + X -
                        24 * (lambda1 * (7 - c1) + lambda2 * (7 - c2) + lambda3 * (7 - c3)))
                eqs.append(term)
            
            # 3 Balance equations
            eqs.append(sum(coeffs_c1[k] * p[k] for k in range(21)) - 4 * c1)
            eqs.append(sum(coeffs_c2[k] * p[k] for k in range(21)) - 4 * c2)
            eqs.append(sum(coeffs_c3[k] * p[k] for k in range(21)) - 4 * c3)
            
            return np.array(eqs)
        
        # Better initial guess based on concentrations
        p_init = np.array([
            c1**4, c2**4, c3**4,  # Pure configs
            4*c1**3*c2, 4*c1**2*c2**2, 2*c1**2*c2**2, 4*c1*c2**3,  # c1-c2 mixed
            4*c1**3*c3, 4*c1**2*c3**2, 2*c1**2*c3**2, 4*c1*c3**3,  # c1-c3 mixed
            4*c2**3*c3, 4*c2**2*c3**2, 2*c2**2*c3**2, 4*c2*c3**3,  # c2-c3 mixed
            4*c1**2*c2*c3, 8*c1*c2*c3**2, 4*c1*c2**2*c3,  # All three
            8*c1*c2**2*c3, 4*c1*c2*c3**2, 8*c1*c2*c3**2
        ])
        p_init = p_init / p_init.sum()
        lambda_init = np.zeros(3)
        initial_guess = np.concatenate([p_init, lambda_init])
        
        # Try optimization with proper method
        from scipy.optimize import least_squares
        
        solution = None
        try:
            # Use least_squares with bounds
            lower_bounds = np.concatenate([np.zeros(21), np.full(3, -100)])
            upper_bounds = np.concatenate([np.ones(21), np.full(3, 100)])
            
            sol = least_squares(
                equations,
                initial_guess,
                bounds=(lower_bounds, upper_bounds),
                method='trf',
                max_nfev=10000,
                verbose=0
            )
            
            if sol.success and np.all(sol.x[:21] > 0):
                solution = sol
        except Exception as e:
            print(f"Optimization failed for T={T}K: {e}")
            continue

        if solution is None or not solution.success or np.any(solution.x[:21] <= 0):
            print(f"Optimization failed for T={T}K.")
            continue
            
        p_values = solution.x[:21]
        
        # Calculate pair probabilities
        pair_probs = {r: {p: 0 for p in req_pairs} for r in [1, 2]}
        
        for r in [1, 2]:
            for (pair, config_k), count in SQUARE_MODEL_N_IJ_K[r].items():
                if pair in pair_probs[r]:
                    pair_probs[r][pair] += count * p_values[config_k - 1]

        # Final Alpha Calculation
        current_temp_alpha[T] = {}
        for name, r in neighbor_map.items():
            current_temp_alpha[T][name] = {}
            for pair, prob in pair_probs[r].items():
                idx1, idx2 = list(pair)
                c_i, c_j = concentrations[idx1 - 1], concentrations[idx2 - 1]
                
                denom = c_i * c_j * SQUARE_MODEL_Z[r] * SQUARE_MODEL_N_STAR
                alpha = 1 - prob / denom if denom != 0 else np.nan
                
                comp_key = components[idx1 - 1] + components[idx2 - 1]
                current_temp_alpha[T][name][f'alpha_{comp_key}'] = alpha

        alpha_values.update(current_temp_alpha)
    
    return alpha_values


def quaternary_sro_pair(fugacity_results, concentrations, pair_names):
    """Calculates SRO for quaternary systems."""
    c1, c2, c3, c4 = concentrations
    alpha_values = {}

    for T, pairs in fugacity_results.items():
        alpha_values[T] = {}
        for neighbor, data in pairs.items():
            try:
                # Ensure all necessary pairs exist
                fugacities = [data[pair_names[i]] for i in range(6)]
            except KeyError as e:
                print(f"Missing fugacity {e} at {T} K")
                continue

            x12, x13, x14, x23, x24, x34 = fugacities

            def equations(vars):
                p = vars
                # p1..p4 are single atom probs? (The logic here seems specific to the model derived)
                # Assuming standard pair cluster expansion logic from original code
                # p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
                eqs = [
                    (p[0] * p[4]) / (p[1]**2) - 0.25 * x12,
                    (p[0] * p[7]) / (p[2]**2) - 0.25 * x13,
                    (p[0] * p[9]) / (p[3]**2) - 0.25 * x14,
                    (p[4] * p[7]) / (p[5]**2) - 0.25 * x23,
                    (p[4] * p[9]) / (p[6]**2) - 0.25 * x24,
                    (p[7] * p[9]) / (p[8]**2) - 0.25 * x34,
                    2 * p[0] + p[1] + p[2] + p[3] - 2 * c1,
                    p[1] + 2 * p[4] + p[5] + p[6] - 2 * c2,
                    p[2] + p[5] + 2 * p[7] + p[6] - 2 * c3,
                    p[3] + p[6] + p[8] + 2 * p[9] - 2 * c3  
                ]
                return eqs

            sol = fsolve(equations, [0.1] * 10)
            if not np.all(np.isfinite(sol)): continue

            # Extract specific probability terms for alpha calc
            # p[1] -> p12, p[2] -> p13, p[3] -> p14, etc.
            alpha_values[T][neighbor] = {
                f'alpha_{pair_names[0]}': 1 - sol[1] / (2 * c1 * c2),
                f'alpha_{pair_names[1]}': 1 - sol[2] / (2 * c1 * c3),
                f'alpha_{pair_names[2]}': 1 - sol[3] / (2 * c1 * c4),
                f'alpha_{pair_names[3]}': 1 - sol[5] / (2 * c2 * c3),
                f'alpha_{pair_names[4]}': 1 - sol[6] / (2 * c2 * c4),
                f'alpha_{pair_names[5]}': 1 - sol[8] / (2 * c3 * c4),
            }
    return alpha_values


def quinary_sro_pair(fugacity_results, concentrations, pair_names):
    """Calculates SRO for quinary systems."""
    c1, c2, c3, c4, c5 = concentrations
    alpha_values = {} 
    
    for T, pairs in fugacity_results.items():
        alpha_values[T] = {}
        for neighbor, data in pairs.items():
            try:
                fugacities = [data[pair_names[i]] for i in range(10)]
            except KeyError as e:
                print(f"Missing fugacity {e} at {T} K")
                continue

            x_vals = fugacities # x12 to x45

            def equations(vars):
                p = vars # 15 variables
                eqs = [
                    (p[0] * p[5]) / (p[1]**2) - 0.25 * x_vals[0],   # x12
                    (p[0] * p[9]) / (p[2]**2) - 0.25 * x_vals[1],   # x13
                    (p[0] * p[12]) / (p[3]**2) - 0.25 * x_vals[2],  # x14
                    (p[0] * p[14]) / (p[4]**2) - 0.25 * x_vals[3],  # x15
                    (p[5] * p[9]) / (p[6]**2) - 0.25 * x_vals[4],   # x23
                    (p[5] * p[12]) / (p[7]**2) - 0.25 * x_vals[5],  # x24
                    (p[5] * p[14]) / (p[8]**2) - 0.25 * x_vals[6],  # x25
                    (p[5] * p[14]) / (p[10]**2) - 0.25 * x_vals[7], # x34 
                    (p[9] * p[12]) / (p[11]**2) - 0.25 * x_vals[8], # x35
                    (p[9] * p[14]) / (p[13]**2) - 0.25 * x_vals[9], # x45
                    # Mass balances
                    2*p[0] + p[1] + p[2] + p[3] + p[4] - 2*c1,
                    p[1] + 2*p[5] + p[6] + p[7] + p[8] - 2*c2,
                    2*p[9] + p[10] + p[11] + p[2] + p[6] - 2*c3,
                    p[10] + 2*p[12] + p[13] + p[3] + p[7] - 2*c3,
                    p[11] + p[13] + 2*p[14] + p[4] + p[8] - 2*c3
                ]
                return eqs

            sol = fsolve(equations, [0.1] * 15)
            if not np.all(np.isfinite(sol)): continue

            alpha_values[T][neighbor] = {
                f'alpha_{pair_names[0]}': 1 - sol[1] / (2 * c1 * c2),
                f'alpha_{pair_names[1]}': 1 - sol[2] / (2 * c1 * c3),
                f'alpha_{pair_names[2]}': 1 - sol[3] / (2 * c1 * c4),
                f'alpha_{pair_names[3]}': 1 - sol[4] / (2 * c1 * c5),
                f'alpha_{pair_names[4]}': 1 - sol[6] / (2 * c2 * c3),
                f'alpha_{pair_names[5]}': 1 - sol[7] / (2 * c2 * c4),
                f'alpha_{pair_names[6]}': 1 - sol[8] / (2 * c2 * c5),
                f'alpha_{pair_names[7]}': 1 - sol[10] / (2 * c3 * c4),
                f'alpha_{pair_names[8]}': 1 - sol[11] / (2 * c3 * c5),
                f'alpha_{pair_names[9]}': 1 - sol[13] / (2 * c4 * c5),
            }
    return alpha_values


def plot_sro_parameters(sro_parameters, alloy_name, output_dir):
    """Generates and saves plots for SRO parameters."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for T, neighbors in sro_parameters.items():
        if not neighbors: continue

        distances = []
        alpha_data = {}

        # Sort neighbor keys (e.g., '1st', '2nd')
        for neighbor in sorted(neighbors.keys()):
            distances.append(neighbor)
            for name, value in neighbors[neighbor].items():
                alpha_data.setdefault(name, []).append(value)

        # Plot setup
        fig, ax = plt.subplots(figsize=(10, 6))
        x_indices = np.arange(len(distances))
        
        for name, values in alpha_data.items():
            ax.plot(x_indices, values, marker='o', label=name)

        ax.set_title(f'{alloy_name} SRO Parameters at {T} K')
        ax.set_xlabel('Nearest Neighbor Distance')
        ax.set_ylabel('SRO Parameter (alpha)')
        ax.set_xticks(x_indices)
        ax.set_xticklabels(distances)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.set_ylim(-0.5, 0.5)
        ax.legend()
        
        filename = output_path / f"{alloy_name}_sro_{T}K.png"
        fig.savefig(filename, dpi=300)
        plt.close(fig)
        print(f"Plot saved: {filename}")


def calculate_sro_parameters(epi_file, sro_output_file, concentrations, number_of_components, fugacity_output_file, nn_list, components, alloy_name, temperatures, alloy_path):
    """
    Backward compatibility wrapper for main_screening.py.
    """
    # 1. Parse EPI Data
    results_dict, parsed_pair_names = parse_epi_file(epi_file, nn_list)
    
    # 2. Calculate Fugacities
    fugacity_results = calculate_fugacity(results_dict, temperatures)
    
    # Save Fugacity
    try:
        with open(fugacity_output_file, 'w') as file:
            json.dump(fugacity_results, file, indent=4)
    except IOError as e:
        print(f"Error writing fugacity file: {e}")

    # 3. Calculate SRO based on component count
    sro_parameters = {}

    prob_output_file = Path(sro_output_file).parent / f'{alloy_name}_pair_probabilities.json'
    
    if number_of_components == 2:
        sro_parameters = binary_sro_pair(concentrations, fugacity_results, parsed_pair_names)
        
    elif number_of_components == 3:
        # The original script had logic to switch between square and pair based on nn_list length
        if len(nn_list) == 2:
            sro_parameters = ternary_sro_square(fugacity_results, concentrations, components, temperatures)
        else:
            sro_parameters = ternary_sro_pair(fugacity_results, concentrations, components, probabilities_output_file=prob_output_file)
            
    elif number_of_components == 4:
        sro_parameters = quaternary_sro_pair(fugacity_results, concentrations, parsed_pair_names)
        
    elif number_of_components == 5:
        sro_parameters = quinary_sro_pair(fugacity_results, concentrations, parsed_pair_names)
    
    else:
        print(f"Error: No model for {number_of_components} components.")
        return {}

    # 4. Save Results
    try:
        with open(sro_output_file, 'w') as file:
            json.dump(sro_parameters, file, indent=4)
            print(f"SRO parameters saved to {sro_output_file}")
    except IOError as e:
        print(f"Error writing SRO file: {e}")

    #plot_sro_parameters(sro_parameters, alloy_name, alloy_path)

    return sro_parameters

