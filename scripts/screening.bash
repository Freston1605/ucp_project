#!/bin/bash
# Script : screening.bash
# 
# Script Description:
# This script automates the molecular screening process using the Autodock VINA tool.
# It is designed to simulate the docking of a multitude of ligands against a specific protein receptor using one 
# configuration file containing the grid parameters.
# It takes three positional arguments in the following order: protein, ligands, and configuration files.

# Arguments:
#   protein: The basename of a protein model file in PDBQT format.
#   ligands: The name of a folder containing all ligands in PDBQT format.
#   configuration (optional): The basename of a configuration file for VINA containing the grid and other docking parameters in txt format.
#       It expects arguments other than --receptor and --ligand but will not enforce it 
#       Unexpected behaviour can occurr as the configuration file mey override the call to vina from this script.

# The script iterates over the ligands folder, finding each appropriate ligand using a wildcard *.pdbqt extension.

# For each ligand, the script invokes VINA with the following arguments:
#   vina --receptor [protein] --config [configuration] --ligand [ligand] --out [output]

# The output is stored in a folder with the following structure: data/output/vina/[protein]/[ligand]/[current_date].

# Additionally, it performs checks to ensure the existence of required directories and files before proceeding, facilitating debugging.

# Usage:
#   bash screening.bash [protein] [ligands] [configuration]

# Example:
#   bash screening.bash 8g8wa flavonoids 
#   bash screening.bash 8j1na flavonoids

# Common folders

# Folder containing all proteins in PDBQT format
PROTEIN_BASE=data/input/pdbqt

# Folder containing all ligands in PDBQT format
LIGAND_BASE=data/output/pdbqt

# Folder containing all configuration files for VINA
CONF_BASE=data/input/conf/vina

# Full path for the protein combining
PROTEIN_PATH=$PROTEIN_BASE/${1}.pdbqt

# Path for the ligand folder containing all ligands
LIGAND_FOLDER=$LIGAND_BASE/${2}

# Full path to list all ligands inside the LIGAND FOLDER
LIGAND_PATH=$LIGAND_FOLDER/*.pdbqt

# Full path for the configuration file
# Check if a third positional argument is supplied
if [ -n "$3" ]; then
    CONF_PATH=$CONF_BASE/${3}.txt
else
    CONF_PATH=$CONF_BASE/${1}.txt
fi

# Grab the current date to label the output folder.
CURRENT_DATE=$(date +%d-%m-%y)

# Define the base directory where output will be stored.
OUTPUT_BASE="/data/output/vina"

# Concatenate the full path for the output folder.
OUTPUT_PATH="$OUTPUT_BASE/${1}/${2}/$CURRENT_DATE"

# Debugging
# Call this script as:
# bash screening.bash [protein] [ligands] > test.txt

# Check if all the PATHS are valid
echo -e "Working with the following paths: \n" \
    "PROTEIN_PATH: $PROTEIN_PATH \n" \
    "LIGAND_PATH: $LIGAND_PATH \n" \
    "CONF_PATH: $CONF_PATH \n" \
    "OUTPUT_PATH: $OUTPUT_PATH \n"

# Check if each path exists and print "OK" if it does
if [ -e "$PROTEIN_PATH" ]; then
    echo "PROTEIN_PATH exists: OK"
else
    echo "PROTEIN_PATH does not exist"
fi

if [ -e "$LIGAND_FOLDER" ]; then
    echo "LIGAND_FOLDER exists: OK"
else
    echo "LIGAND_FOLDER does not exist"
fi

if [ -e "$CONF_PATH" ]; then
    echo "CONF_PATH exists: OK"
else
    echo "CONF_PATH does not exist"
fi

if [ -e "$OUTPUT_PATH" ]; then
    echo "OUTPUT_PATH exist"
else
    echo "OUTPUT_PATH exists already: OK"
fi



# Loop over the files with the '.pdbqt' extension in the LIGAND_PATH folder
for file in $LIGAND_PATH; do
	# Check if the file exists
    if [ -f "$file" ]; then
    # Debugging
    # Print to console the ligand processed and the call to VINA
    base_name=$(basename $file .pdbqt)
    echo "Processing ligand: $file"
    echo "mkdir -p $OUTPUT_PATH/$base_name"
    echo "vina \\
    --receptor $PROTEIN_PATH \\
    --config $CONF_PATH \\
    --ligand $file \\
    --out $OUTPUT_PATH/$base_name/$base_name.out"
	# Molecular Docking
    # mkdir -p $OUTPUT_PATH/$base_name
	# vina\
	# --receptor $PROTEIN_PATH \
	# --config $CONF_PATH \
	# --ligand $file \
	# --out $OUTPUT_PATH/$base_name/$base_name/$base_name.out
    fi
done
