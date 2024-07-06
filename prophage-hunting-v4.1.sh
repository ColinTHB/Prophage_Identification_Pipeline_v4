#!/bin/bash
set -euo pipefail

# Initialize Conda
source ~/miniconda3/etc/profile.d/conda.sh

# Function to display help message
show_help() {
    echo "Usage: $0 -i input_fasta_file"
    echo "  -i input_fasta_file   Specify the input genome fasta sequence file."
}

# Parse command-line arguments
while getopts ":i:h" opt; do
    case ${opt} in
        i )
            input_fasta=${OPTARG}
            ;;
        h )
            show_help
            exit 0
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            show_help
            exit 1
            ;;
        : )
            echo "Invalid option: -$OPTARG requires an argument" 1>&2
            show_help
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

# Check if the input fasta file is provided
if [ -z "$input_fasta" ]; then
    echo "Error: Input fasta file is required."
    show_help
    exit 1
fi

# Define the base name for the output directories and files
base_name=$(basename "$input_fasta" .fasta)

# Step 1: Run VirSorter with the input fasta file
echo "Running VirSorter..."
conda activate virsorter2_env
virsorter run -w ${base_name}_virsorter_output -i "$input_fasta" --include-groups "dsDNAphage,ssDNA,RNA" -j 10 -d /mnt/colinb/databases/virsorter2/
if [ $? -ne 0 ]; then
    echo "VirSorter failed."
    conda deactivate
    exit 1
fi
conda deactivate

# Step 1.1: Rewrite the first column of the final-viral-boundary.tsv file with first column of the final-viral-score.tsv issue with formatting of former with names

# Extract the first column from final-viral-score.tsv
cut -f1 ${base_name}_virsorter_output/final-viral-score.tsv > ${base_name}_virsorter_output/temp-column.txt

# Combine the extracted column with the rest of final-viral-boundary.tsv
paste ${base_name}_virsorter_output/temp-column.txt <(cut -f2- ${base_name}_virsorter_output/final-viral-boundary.tsv) > ${base_name}_virsorter_output/final-viral-boundary-new.tsv

# Clean up temporary file
rm ${base_name}_virsorter_output/temp-column.txt

# Step 2: Select prophage regions of interest
echo "Selecting strong prophage regions..."
awk '($2 >= 0.85 && $8 > 2 && $6 ~ /dsDNAphage/)' ${base_name}_virsorter_output/final-viral-score.tsv | cut -f 1 > ${base_name}_virsorter_output/virsorter_strong.ids
if [ $? -ne 0 ]; then
    echo "Selection of strong prophage regions failed."
    exit 1
fi

# Step 3: Grab prophage final-viral-score.tsv details
echo "Grabbing final viral score details..."
grep -Ff ${base_name}_virsorter_output/virsorter_strong.ids ${base_name}_virsorter_output/final-viral-score.tsv > ${base_name}_virsorter_output/selected-final-score.tsv
if [ $? -ne 0 ]; then
    echo "Grabbing final viral score details failed."
    exit 1
fi

# Step 4: Grab prophage genomes of interest
echo "Extracting prophage genomes of interest..."
conda activate pullseq_env
pullseq -i ${base_name}_virsorter_output/final-viral-combined.fa -n ${base_name}_virsorter_output/virsorter_strong.ids > ${base_name}_virsorter_output/selected_prophage_sequences.fa
if [ $? -ne 0 ]; then
    echo "Extraction of prophage genomes failed."
    conda deactivate
    exit 1
fi
conda deactivate

# Step 5: Annotate prophage genomes with pharokka
echo "Annotating prophage genomes with pharokka..."
conda activate pharokka_env
pharokka.py -i ${base_name}_virsorter_output/selected_prophage_sequences.fa -o ${base_name}_virsorter_output/pharokka_results_2024 -d /mnt/colinb/databases/pharokka/ -t 10
if [ $? -ne 0 ]; then
    echo "Pharokka annotation failed."
    conda deactivate
    exit 1
fi
conda deactivate

# Step 6: Extract prophage gene/protein details
echo "Extracting prophage gene/protein details..."
grep -E "portal protein|terminase large subunit|major head proteins" ${base_name}_virsorter_output/pharokka_results_2024/pharokka.gff | awk -F'\t' '{
    sequence_name = $1
    match($9, /ID=([^;]+)/, id_match)
    id = id_match[1]
    start = $4
    end = $5
    match($9, /product=([^;]+)/, product_match)
    product = product_match[1]
    print sequence_name "\t" id "\t" start "\t" end "\t" product
}' > ${base_name}_virsorter_output/improved_output.tsv
if [ $? -ne 0 ]; then
    echo "Extraction of prophage gene/protein details failed."
    exit 1
fi

# Step 7: Merge selected-final-score.tsv with improved_output.tsv
echo "Merging selected final score with improved output..."
awk 'BEGIN {
    FS = "\t"; OFS = "\t"
    while (getline < "'"${base_name}_virsorter_output/selected-final-score.tsv"'") {
        selected_map[$1] = $0
    }
}
{
    split($0, arr, "\t")
    seqname = arr[1]
    if (seqname in selected_map) {
        print $0, selected_map[seqname]
    } else {
        print $0
    }
}' ${base_name}_virsorter_output/improved_output.tsv > ${base_name}_virsorter_output/merged-output.tsv
if [ $? -ne 0 ]; then
    echo "Merging final score with improved output failed."
    exit 1
fi

# Step 8: Add trim bp start and end from final-viral-boundary.tsv to merged-output.tsv 
echo "Adding trim bp start and end to merged output..." 
boundary_file="${base_name}_virsorter_output/final-viral-boundary-new.tsv" 
merge_output_file="${base_name}_virsorter_output/merged-output.tsv" 
final_output_file="${base_name}_virsorter_output/final-output.tsv"

if [ ! -f "$boundary_file" ]; then 
    echo "Error: $boundary_file does not exist." 
    exit 1 
fi 
 
if [ ! -f "$merge_output_file" ]; then 
    echo "Error: $merge_output_file does not exist." 
    exit 1 
fi 

awk 'BEGIN { 
    FS = "\t"; OFS = "\t"; 
    
    # Read boundary file and store trim_bp_start and trim_bp_end in arrays
    while (getline < "'"$boundary_file"'") { 
        seqname = $1; 
        trim_bp_start[seqname] = $4; 
        trim_bp_end[seqname] = $5; 
    } 
    close("'"$boundary_file"'"); 
} 
{
    # Match the sequence name pattern in the current line of merged output file
    match($0, /NODE_[0-9]+_length_[0-9]+_cov_[0-9]+\|\|full|NODE_[0-9]+_length_[0-9]+_cov_[0-9]+\|\|[0-9]+_partial/); 
    seqname = substr($0, RSTART, RLENGTH); 
    
    # Print the line with trim_bp_start and trim_bp_end or "NA" if not found
    if (seqname in trim_bp_start && seqname in trim_bp_end) { 
        print $0, trim_bp_start[seqname], trim_bp_end[seqname]; 
    } else { 
        print $0, "NA", "NA"; 
    } 
}' "$merge_output_file" > "$final_output_file"

if [ $? -ne 0 ]; then 
    echo "Adding trim bp start and end failed." 
    exit 1 
fi 

echo "Trim bp start and end added successfully to $final_output_file."

# Step 9: Add headers to final-output.tsv
echo "Adding headers to final output..."
echo -e "seq_name\tlocus_tag\tstart\tend\tproduct\tseq_name\tdsDNAphage\tssDNA\tRNA\tmax_score\tmax_score_group\tlength\thallmark\tviral\tcellular\ttrim_bp_start\ttrim_bp_end" | cat - ${base_name}_virsorter_output/final-output.tsv > ${base_name}_virsorter_output/temp && mv ${base_name}_virsorter_output/temp ${base_name}_virsorter_output/final-output.tsv
if [ $? -ne 0 ]; then
    echo "Adding headers failed."
    exit 1
fi

echo "Pipeline completed successfully."
