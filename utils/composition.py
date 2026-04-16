import re
import pandas as pd


# Interactive input for alloy composition
def get_composition_from_user(default=False):
    if default:
        return {"Fe": 30, "Ni": 30, "Cr": 40}
    else:
        print("\nEnter the alloy composition (e.g., Fe 30, Ni 30, Cr 40):")
    composition_input = input("Composition: ")
    try:
        # Convert input string into a dictionary
        comp = {item.split()[0]: int(item.split()[1]) for item in composition_input.split(",")}
        return comp
    except Exception as e:
        print(f"Error parsing composition: {e}")
        return {}

# Function to parse composition strings from Excel
def parse_composition_string(comp_string):
    element_pattern = re.compile(r'([A-Z][a-z]?)(\d+)')
    matches = element_pattern.findall(comp_string)
    return {element: int(count) for element, count in matches}

# Load compositions from Excel
def load_compositions_from_excel():
    path = input("Enter the path to the .xls file: ")
    df = pd.read_excel(path, sheet_name="Alloy systems")
    compositions = df.iloc[0:, 0].tolist() # Skip the header
    return [parse_composition_string(comp) for comp in compositions if isinstance(comp, str)]

# Calculate the number of components in alloy
def number_of_components(composition):
    n = len(composition.keys())
    print(f"Number of components in alloy: {n}")
    return n

# Calculate the number of pairs in alloy
def number_of_pairs(composition):
    n = number_of_components(composition)
    # Formula for number of unique pairs: n(n-1)/2
    pairs = int((n * (n - 1)) / 2)
    print(f"Number of unique pairs in alloy: {pairs}")
    return pairs
def get_alloy_name_and_concentrations(composition):
    # Create alloy name by concatenating elements and their atomic percentages
    alloy_name = ''.join([f"{element}{int(concentration)}" for element, concentration in composition.items()])
    
    # Convert percentages to decimal concentrations
    concentrations = [concentration/100 for concentration in composition.values()]
    
    print(f"Alloy name: {alloy_name}")
    print(f"Concentrations: {concentrations}")
    
    return alloy_name, concentrations
           