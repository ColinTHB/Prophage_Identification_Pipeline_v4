import os
import argparse

def split_gff(input_file, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    current_sequence = None
    current_lines = []
    sequence_count = 1

    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('##sequence-region'):
                # Skip these lines in GFF3 format, as they are not feature records
                continue
            elif line.startswith('#'):
                # Skip comment lines
                continue
            elif line == '':
                # Skip empty lines
                continue
            else:
                # Split the line into columns based on tab delimiter
                columns = line.split('\t')
                if len(columns) < 9:
                    continue  # Skip malformed lines
                
                # Extract sequence identifier from the first column
                sequence_id = columns[0]

                # Check if we encountered a new sequence identifier
                if sequence_id != current_sequence:
                    # Write out the previous sequence set if it exists
                    if current_sequence is not None:
                        output_file = os.path.join(output_dir, f"pharokka_{sequence_count}.gff")
                        with open(output_file, 'w') as outfile:
                            outfile.writelines(current_lines)
                            print(f"Written {output_file}")
                        sequence_count += 1

                    # Reset for the new sequence identifier
                    current_sequence = sequence_id
                    current_lines = []

                # Accumulate lines for the current sequence identifier
                current_lines.append(line + '\n')

        # Write out the last sequence set
        if current_sequence is not None:
            output_file = os.path.join(output_dir, f"pharokka_{sequence_count}.gff")
            with open(output_file, 'w') as outfile:
                outfile.writelines(current_lines)
                print(f"Written {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a multi-sequence GFF file into individual files named pharokka_1.gff, pharokka_2.gff, etc.")
    parser.add_argument("input_file", type=str, help="Path to the input multi-sequence GFF file.")
    parser.add_argument("output_dir", type=str, help="Directory where the output files will be saved.")

    args = parser.parse_args()
    
    split_gff(args.input_file, args.output_dir)
