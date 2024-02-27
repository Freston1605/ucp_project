import os
import subprocess
from datetime import datetime
import sys
import glob
import argparse
import logging


def print_working_paths(
    protein_path, ligand_dir, conf_path, out_path, ligand_base, ligand_full
):
    """
    Print working paths and variables.

    Parameters:
    - protein_path (str): Path to the protein receptor file.
    - ligand_path (str): Path to the directory containing ligand files.
    - conf_path (str): Path to the Vina configuration file.
    - out_path (str): Output path for the results.
    - ligand_base (str): Base name for ligand files (wildcard pattern).
    - ligand_full (str): Full path to all ligand files.
    """
    logging.info("Working with the following paths:")
    logging.info(f"PROTEIN_PATH: {protein_path}")
    logging.info(f"LIGAND_PATH: {ligand_dir}")
    logging.info(f"CONF_PATH: {conf_path}")
    logging.info(f"OUT_PATH: {out_path}")
    logging.info(f"LIGAND_BASE: {ligand_base}")
    logging.info(f"LIGAND_FULL: {ligand_full}\n")


def process_ligand(protein_path, out_path, file):
    """
    Process a ligand file.

    Parameters:
    - protein_path (str): Path to the protein receptor file.
    - out_path (str): Output path for the results.
    - file (str): Full path to the ligand file.
    """
    base_name = os.path.basename(file).split(".pdbqt")[0]
    logging.info(f"Processing ligand: {base_name}")

    # Create output directory
    out_dir = os.path.join(out_path, base_name)
    os.makedirs(out_dir, exist_ok=True)

    # Run Vina
    subprocess.run(
        [
            "vina",
            "--receptor",
            protein_path,
            "--config",
            conf_path,
            "--ligand",
            file,
            "--out",
            os.path.join(out_dir, f"{base_name}.out"),
        ]
    )


if __name__ == "__main__":
    """
    Molecular Screening Script

    This script automates molecular screening using the Vina program in a SLURM environment. It processes ligands
    against a specified protein receptor, creating output directories for each ligand with the results.

    Usage:
    python script_name.py protein_name ligands_directory

    Example:
    python molecular_screening.py t2r4 stevia

    This example runs molecular screening for the protein 't2r4' using ligands from the 'stevia' directory.
    """
    # Argument parsing
    parser = argparse.ArgumentParser(description="Molecular screening script.")
    parser.add_argument("protein", help="Name of the protein")
    parser.add_argument("ligands", help="Name of the ligands directory")
    args = parser.parse_args()

    # Base directories
    base_work = "/work/waldo.acevedo/sg_block/screening"
    base_scratch = "/scratch/waldo.acevedo/sg_block"

    # Current date
    current_date = datetime.now().strftime("%d-%m-%y")

    # Directories and files
    protein_path = f"{base_work}/alphafold/{args.protein}/{args.protein}.pdbqt"
    ligand_path = f"{base_work}/ligands/{args.ligands}"
    conf_path = f"{base_work}/alphafold/{args.protein}/conf_{args.protein}_stevia.txt"
    out_path = f"{base_scratch}/{args.ligands}/{args.protein}/{current_date}"
    
    ligand_pattern = "ligand_*.pdbqt"
    ligand_full = os.path.join(ligand_path, ligand_pattern)

    # Set up logging
    logging.basicConfig(filename="screening_log.txt", level=logging.INFO)

    # Print paths
    print_working_paths(
        protein_path, ligand_path, conf_path, out_path, ligand_pattern, ligand_full
    )

    # Loop through ligands
    for ligand_file in glob.glob(ligand_full):
        process_ligand(protein_path, out_path, ligand_file)
