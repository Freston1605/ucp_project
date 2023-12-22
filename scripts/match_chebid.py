from rdkit import Chem
from rdkit.Chem import AllChem
import csv
import os

def read_chebi_ids_from_csv(csv_file_path):
    """Read ChEBI IDs from a CSV file.

    Args:
        csv_file_path (str): Path to the CSV file.

    Returns:
        list: List of ChEBI IDs.
    """
    chebi_ids = []
    with open(csv_file_path, "r", encoding="utf-8") as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader)  # Skip the header row
        chebi_ids = [row[0] for row in csv_reader]
    return chebi_ids

def read_sdf_file(input_sdf_file):
    """Read an SDF file using RDKit.

    Args:
        input_sdf_file (str): Path to the input SDF file.

    Returns:
        Chem.SDMolSupplier: RDKit molecule supplier.
    """
    return Chem.SDMolSupplier(input_sdf_file)

def write_matched_molecules_to_sdf(output_sdf_file, matched_molecules):
    """Write matched molecules to an SDF file.

    Args:
        output_sdf_file (str): Path to the output SDF file.
        matched_molecules (list): List of RDKit molecules containing matched molecules.
    """
    # Use a temporary molecule supplier to create the output SDF file
    with Chem.SDWriter(output_sdf_file) as writer:
        for mol in matched_molecules:
            writer.write(mol)

def main(input_sdf, input_csv, output="matched_molecules.sdf"):
    """Main function to process molecules and write matched molecules to an SDF file.

    Args:
        input_sdf (str): Path to the input SDF file.
        input_csv (str): Path to the input CSV file containing ChEBI IDs.
        output (str, optional): Path to the output SDF file. Defaults to "matched_molecules.sdf".
    """

    # Read ChEBI IDs from the CSV file
    chebi_ids = read_chebi_ids_from_csv(input_csv)

    # Read original SDF file
    mol_supplier = read_sdf_file(input_sdf)

    # Create a new SDF file to store matching molecules
    output = os.path.join(output)

    # Check and create the output SDF files if it does not exist
    if not os.path.exists(output):
        with open(output, "w", encoding="utf-8") as file:
            file.write("Flavonoids" + "\n")
            print("Output file created:", output)

    # Write matched molecules to the new SDF file
    matched_molecules = [mol for mol in mol_supplier if (mol is not None and mol.GetProp("ChEBI ID").split(":")[1] in chebi_ids)]
    write_matched_molecules_to_sdf(output, matched_molecules)

    print(
        f"Finished processing molecules. Matched molecules written to '{output}'."
    )

# Execution


# Main paths
INPUT_PATH = "data/input"
OUTPUT_PATH = "data/output"

# CSV containing ChEBI IDs
INPUT_CSV = "csv/chebi_ids.csv"
INPUT_CSV = os.path.join(OUTPUT_PATH, INPUT_CSV)

# SDF containing al molecules
INPUT_SDF = "sdf/ChEBI_complete.sdf"
INPUT_SDF = os.path.join(INPUT_PATH, INPUT_SDF)

# OUTPUT SDF containing all matching molecules
OUTPUT_SDF = os.path.join(OUTPUT_PATH, "sdf/matched_molecules.sdf")

main(INPUT_SDF, INPUT_CSV, OUTPUT_SDF)
