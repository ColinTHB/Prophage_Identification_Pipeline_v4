#!/bin/bash

# Check if the user provided the required argument
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_virsorter_output>"
    exit 1
fi

# Warn the user if BASE_DIR ends with a trailing slash
if [[ "$1" == */ ]]; then
    echo "Warning: The provided path should not end with a trailing slash (/)."
    echo "This may affect proper file operations within the script."
    exit 1
fi

# Initialize Conda (assuming this is correct and intentional)
source ~/miniconda3/etc/profile.d/conda.sh

# Activate Conda environment with bedtools and biopython
conda activate bedtools_env

# Define the base directory based on the provided argument
BASE_DIR="$1"

# Split the pharokka gff file into its subfiles
python ./processing_scripts_v3/split_gff.py "$BASE_DIR/pharokka_results_2024/pharokka.gff" "$BASE_DIR/pharokka_results_2024/"

# Define the output directory
output_directory="$BASE_DIR/pharokka_results_2024/"

# Process each gff file in the target directory
for input_file in "$BASE_DIR/pharokka_results_2024/pharokka_"*.gff; do
    if [[ "$input_file" != "$BASE_DIR/pharokka_results_2024/pharokka.gff" ]] && [ -f "$input_file" ]; then
        python ./processing_scripts_v3/process_pharokka_gff_for_prophage_boundaries-v3.py "$input_file" "$output_directory"
    fi
done

# Call script to extract prophage boundaries and write to file - dynamically process all gff output files
for input_file in "$BASE_DIR/pharokka_results_2024/pharokka_"*_output.txt; do
    if [[ "$input_file" == *"_output.txt" ]] && [ -f "$input_file" ]; then
        python ./processing_scripts_v3/extract_prophage_boundaries-v5.py "$input_file"
    fi
done

# Remove unwanted files
rm "$BASE_DIR/pharokka_results_2024/pharokka_minced_output.txt"
rm "$BASE_DIR/pharokka_results_2024/pharokka_aragorn_output.txt"

# Split the pharokka gbk file into its subfiles
python ./processing_scripts_v3/split_genbank.py "$BASE_DIR/pharokka_results_2024/pharokka.gbk" "$BASE_DIR/pharokka_results_2024/"

# Call script to extract prophage and write to directory where the executing script resides
for input_file in "$BASE_DIR/pharokka_results_2024/pharokka_"*.gbk; do
    output_file="${input_file%.gbk}_output_prophage_boundaries.tsv"
    if [ -f "$input_file" ]; then
        ./processing_scripts_v3/subset_genbank-v2.sh "$input_file" "$output_file"
    fi
done

# Remove unwanted files
rm "$BASE_DIR/pharokka_results_2024/pharokka_minced_output_prophage_boundaries.tsv"
rm "$BASE_DIR/pharokka_results_2024/pharokka_aragorn_output_prophage_boundaries.tsv"

# Concatenate prophage_boundaries.txt 
awk 'FNR==1{fname=FILENAME; sub(".*/", "", fname)} {print fname, $0}' OFS='\t' "$BASE_DIR/pharokka_results_2024/pharokka_"*_output_prophage_boundaries.tsv > ./prophage_boundaries_combined.tsv

# Add prophage details to final-output.tsv
./processing_scripts_v3/add-prophage-boundaries-2-tsv-v2.sh -1 "$BASE_DIR/final-output.tsv" -2 ./prophage_boundaries_combined.tsv -o "${BASE_DIR}_final-output-with-prophage.tsv"

# Copy pharokka_*_subset.gbk files to the current directory and rename them with $BASE_DIR prefix
for gbk_file in "$BASE_DIR/pharokka_results_2024/pharokka_"*_subset.gbk; do
    if [ -f "$gbk_file" ]; then
        filename=$(basename "$gbk_file")
        cp "$gbk_file" "./${BASE_DIR}_${filename}"
    fi
done

# Remove unwanted files
rm ./prophage_boundaries_combined.tsv


