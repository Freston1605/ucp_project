from rdkit import Chem
from rdkit.Chem import SDMolSupplier
from pathlib import Path


def separate_molecules(sdf_file_path=str, output_dir=str):
    """Separate molecules from an SDF file and save each molecule as a separate SDF file.
    Every molecule is named after its ChEBI Name and, if None, after its index in the loop.

    Args:
        sdf_file_path (str): Path to the input SDF file.
        output_dir (str): Path to the output directory where separate SDF files will be saved.
    """
    # Create the output directory if it doesn't exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Use SDMolSupplier to read the SDF file
    suppl = SDMolSupplier(sdf_file_path)

    # Iterate over each molecule in the SDF file
    for index, mol in enumerate(suppl):
        if mol is not None:
            # Extract the ChEBI ID (assuming it's stored in the 'ChEBI ID' field)
            chebi_id = mol.GetProp('ChEBI ID')

            # Split and extract the ChEBI ID from the field value
            chebi_id = chebi_id.split(':')[1] if chebi_id else None
            
            # If the ChEBI ID is not available, generate a default name
            if not chebi_id:
                chebi_id = f'molecule_{index + 1}'

            # Generate a unique filename for each molecule
            output_file_path = output_dir / f'{chebi_id}.sdf'

            # Write the molecule to a separate SDF file
            writer = Chem.SDWriter(str(output_file_path))
            writer.write(mol)
            writer.close()

            print(f'Molecule "{chebi_id}" saved to {output_file_path}')

# Execution

sdf_file_path = 'data/output/sdf/matched_molecules.sdf'
output_directory = 'data/output/sdf/flavonoids'

separate_molecules(sdf_file_path, output_directory)