from rdkit import Chem
import csv
import glob
from os import path


# Listing all SDF files on the folder
directory_path = "data/output/sdf/flavonoids"  # Folder containing the files
file_pattern = "*.sdf"  # File extension
file_list = glob.glob(path.join(directory_path, file_pattern))

# Empty list for all molecules
data = []

# Loop trough the files in the folder
for file in file_list:
    mol_supplier = Chem.SDMolSupplier(file)
    mol = next(mol_supplier)
    if mol is not None:
        # Register important molecular data for reference
        chebi_id = mol.GetProp("ChEBI ID") if mol.HasProp("ChEBI ID") else None
        chebi_id = chebi_id.split(":")[1] if chebi_id else None
        chebi_name = mol.GetProp("ChEBI Name") if mol.HasProp("ChEBI Name") else None
        smiles = mol.GetProp("SMILES") if mol.HasProp("SMILES") else None
        pubchem_id = (
            mol.GetProp("PubChem Database Links")
            if mol.HasProp("PubChem Database Links")
            else None
        )  # Write a dictionary with the important data
        mol_dict = {
            "ChEBI ID": chebi_id,
            "ChEBI Name": chebi_name,
            "SMILES": smiles,
            "PubChem SID": pubchem_id,
        }
        # Append the dictionary to the list
        data.append(mol_dict)
        print(f"Appended to the list the molecule: {mol_dict}")

# Output CSV file
csv_file_path = "data/output/csv/flavonoids_data.csv"

# Write the data list to the specified CSV
with open(csv_file_path, "w", encoding="utf-8", newline="") as csv_file:
    field_names = ["ChEBI ID", "ChEBI Name", "SMILES", "PubChem SID"]
    writer = csv.DictWriter(csv_file, fieldnames=field_names)

    writer.writeheader()  # Write the header
    writer.writerows(data)  # Write the data

print(f"CSV file created at: {csv_file_path}")
