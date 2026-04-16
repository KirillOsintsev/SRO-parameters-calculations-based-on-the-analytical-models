#!/usr/bin/env python3
# import os
import json
import itertools
from pathlib import Path

def get_component_pairs(components):
    """Generate unique pairs of components."""
    return list(itertools.combinations(sorted(components), 2))

def get_user_input_nearest_neighbors():
    """Ask user how many nearest neighbors to calculate (1-4)."""
    while True:
        try:
            n = int(input("How many nearest neighbors to calculate? (1-4): "))
            if 1 <= n <= 4:
                return [f"{i}st" if i == 1 else f"{i}nd" if i == 2 else f"{i}rd" if i == 3 else f"{i}th" 
                       for i in range(1, n + 1)]
            print("Please enter a number between 1 and 4")
        except ValueError:
            print("Please enter a valid number")

def get_lattice_constant(components, database_path="../data/input/database.json"):
    """Get lattice constant either from database average or user input."""
    while True:
        choice = input("Calculate lattice constant from database (d) or provide custom value (c)? [d/c]: ").lower()
        
        if choice == 'd':
            try:
                with open(database_path, 'r') as f:
                    db = json.load(f)
                constants = [db[component]["lattice_constant"] for component in components]
                return sum(constants) / len(constants)
            except (FileNotFoundError, KeyError) as e:
                print(f"Error accessing database: {e}")
                print("Falling back to custom input...")
            
        if choice == 'c':
            try:
                value = float(input("Enter lattice constant value: "))
                return value
            except ValueError:
                print("Please enter a valid number")

def get_lattice_constant_auto():
                value = float(input("Enter lattice constant value: "))
                return value               
                

def get_ncells():
    """Get number of unit cells from user or use default."""
    while True:
        choice = input("Use default number of unit cells (5)? [y/n]: ").lower()
        if choice == 'y':
            return 5
        elif choice == 'n':
            try:
                value = int(input("Enter number of unit cells: "))
                return value
            except ValueError:
                print("Please enter a valid number")

def get_phase():
    """Get crystal phase from user."""
    while True:
        phase = input("Enter crystal phase (fcc/bcc): ").lower()
        if phase in ['fcc', 'bcc']:
            return phase
        print("Please enter either 'fcc' or 'bcc'")

def create_folder_structure(alloy_path, components, nn_list):
    """Create the folder structure for calculations."""
    base_path = Path(alloy_path)
    base_path.mkdir(parents=True,exist_ok=True)
    
    pairs = get_component_pairs(components)
    
    for pair in pairs:
        pair_name = ''.join(pair)
        pair_path = base_path / pair_name
        pair_path.mkdir(exist_ok=True)
        
        for nn in nn_list:
            nn_path = pair_path / nn
            nn_path.mkdir(exist_ok=True)
            
    
    return pairs

