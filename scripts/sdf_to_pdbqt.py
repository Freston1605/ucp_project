import os
from openbabel import pybel

def convert_sdf_to_pdbqt(input_folder, output_folder):
    """Converts all of the SD files in a folder and outputs them in pdbqt format in another folder.

    Args:
        input_folder (str): The path of the input folder containing the SD files to be converted.
        output_folder (str): The path of the output folder that will contain the pdbqt files.
    """
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Iterate through the SDF files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith(".sdf"):
            input_path = os.path.join(input_folder, filename)
    
            # Generate the corresponding output path with .pdbqt extension
            output_filename = os.path.splitext(filename)[0] + ".pdbqt"
            output_path = os.path.join(output_folder, output_filename)
            
            # Open the SDF file
            mol_supplier = pybel.readfile("sdf", input_path)
            
            # Create a PDBQT file
            for mol in mol_supplier:
                    mol.write(format="pdbqt", filename=output_path, overwrite=True)
                    

INPUT_FOLDER = "data/output/sdf/flavonoids"
OUTPUT_FOLDER = "data/output/pdbqt/flavonoids"
convert_sdf_to_pdbqt(INPUT_FOLDER, OUTPUT_FOLDER)
