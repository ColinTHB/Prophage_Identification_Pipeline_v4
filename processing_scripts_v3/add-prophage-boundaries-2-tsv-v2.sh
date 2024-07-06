#!/bin/bash

# Default output file
output_file="final-output-with-prophage.tsv"

# Function to display usage
usage() {
    echo "Usage: $0 -1 <final_output.tsv> -2 <prophage_boundaries.tsv> -o <output_file>"
    exit 1
}

# Parse command line options
while getopts ":1:2:o:" opt; do
    case $opt in
        1) final_output="$OPTARG";;
        2) prophage_boundaries="$OPTARG";;
        o) output_file="$OPTARG";;
        \?) echo "Invalid option -$OPTARG" >&2; usage;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage;;
    esac
done

# Check if required options are provided
if [ -z "$final_output" ] || [ -z "$prophage_boundaries" ]; then
    echo "Error: Missing required input file(s)." >&2
    usage
fi

# Check if input files exist
if [ ! -f "$final_output" ]; then
    echo "Error: $final_output not found." >&2
    exit 1
fi

if [ ! -f "$prophage_boundaries" ]; then
    echo "Error: $prophage_boundaries not found." >&2
    exit 1
fi

# Process final-output.tsv
echo -e "eq_name\tlocus_tag\tstart\tend\tproduct\tseq_name\tdsDNAphage\tssDNA\tRNA\tmax_score\tmax_score_group\tlength\thallmark\tviral\tcellular\ttrim_bp_start\ttrim_bp_end\tprophage_left\tprophage_right" > "$output_file"

# Read prophage boundaries into an associative array
declare -A prophage_left prophage_right
while IFS=$'\t' read -r filename seq_name start end; do
    prophage_left["$seq_name"]=$start
    prophage_right["$seq_name"]=$end
done < "$prophage_boundaries"

# Read each line of final-output.tsv
tail -n +2 "$final_output" | while IFS=$'\t' read -r seq_name locus_tag start end product seq_name dsDNAphage ssDNA RNA max_score max_score_group length hallmark viral cellular trim_bp_start trim_bp_end; do
    # Ensure trim_bp_start is numeric
    if ! [[ $trim_bp_start =~ ^[0-9]+$ ]]; then
        echo "Warning: trim_bp_start is not numeric for line $seq_name. Skipping..."
        continue
    fi

    # Retrieve prophage boundaries for current seq_name
    current_prophage_left="${prophage_left[$seq_name]}"
    current_prophage_right="${prophage_right[$seq_name]}"

    if [ -z "$current_prophage_left" ] || [ -z "$current_prophage_right" ]; then
        echo "Warning: Prophage boundaries not found for $seq_name. Skipping..."
        continue
    fi

    # Calculate prophage left and right based on retrieved values
    prophage_left_sum=$((current_prophage_left + trim_bp_start))
    prophage_right_sum=$((current_prophage_right + trim_bp_start))

    # Append the calculated values to the output file
    echo -e "$seq_name\t$locus_tag\t$start\t$end\t$product\t$seq_name\t$dsDNAphage\t$ssDNA\t$RNA\t$max_score\t$max_score_group\t$length\t$hallmark\t$viral\t$cellular\t$trim_bp_start\t$trim_bp_end\t$prophage_left_sum\t$prophage_right_sum" >> "$output_file"
done

echo "Output written to $output_file"
