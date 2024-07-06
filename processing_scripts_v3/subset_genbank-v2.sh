#!/bin/bash

# Check for required arguments
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_gbk> <prophage_boundaries_txt>"
  exit 1
fi

# Input arguments
INPUT_GBK="$1"
PROPHAGE_BOUNDARIES_TXT="$2"

# Determine output directory and file name
INPUT_DIR=$(dirname "$(realpath "$INPUT_GBK")")
OUTPUT_GBK="$INPUT_DIR/$(basename "$INPUT_GBK" .gbk)_subset.gbk"

# Directory of prophage_boundaries.txt
COORDINATES_FILE="$(realpath "$PROPHAGE_BOUNDARIES_TXT")"

if [ ! -f "$COORDINATES_FILE" ]; then
  echo "Coordinates file not found: $COORDINATES_FILE"
  exit 1
fi

# Extract coordinates from the file
START=$(awk '{print $2}' "$COORDINATES_FILE")
END=$(awk '{print $3}' "$COORDINATES_FILE")

if [ -z "$START" ] || [ -z "$END" ]; then
  echo "Failed to extract coordinates from the file"
  exit 1
fi

# Extract the sequence region using bedtools and BioPython
python3 - <<EOF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

input_gbk = "$INPUT_GBK"
start = int("$START")
end = int("$END")
output_gbk = "$OUTPUT_GBK"

# Read the input GenBank file
record = SeqIO.read(input_gbk, "genbank")

# Extract the sequence within the given coordinates
subset_seq = record.seq[start-1:end]

# Create a new record for the subset sequence
new_record = SeqRecord(subset_seq, id=record.id, description=record.description)

# Transfer annotations (if any) and update feature coordinates
new_record.annotations = record.annotations

# Update features
new_features = []
for feature in record.features:
    if start <= feature.location.end and end >= feature.location.start:
        new_start = max(0, feature.location.start - start + 1)
        new_end = min(len(subset_seq), feature.location.end - start + 1)
        new_location = FeatureLocation(new_start, new_end, strand=feature.location.strand)
        new_feature = SeqFeature(location=new_location, type=feature.type, qualifiers=feature.qualifiers)
        new_features.append(new_feature)

new_record.features = new_features

# Write the subset sequence to the output GenBank file
SeqIO.write(new_record, output_gbk, "genbank")
EOF

echo "Subset GenBank file created: $OUTPUT_GBK"
