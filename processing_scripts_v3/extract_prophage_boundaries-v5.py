import re
import sys
import os

def extract_boundaries(text):
    # Regular expression to extract boundary coordinates from the specific statements
    regex = r'Boundaries for ([^:]+): (\d+) to (\d+)'
    matches = re.findall(regex, text)
    return [(identifier, int(start), int(end)) for identifier, start, end in matches]

def select_greatest_range(boundaries):
    if not boundaries:
        return None
    
    # Filter out pairs where the distance is greater than 60,000
    valid_boundaries = [(identifier, start, end) for identifier, start, end in boundaries if abs(end - start) <= 60000]
    
    if not valid_boundaries:
        return None
    
    # Select the pair with the greatest range
    greatest_range_pair = max(valid_boundaries, key=lambda x: x[2] - x[1])
    return greatest_range_pair

def extract_node_details(text):
    # Regular expression to pull NODE details based on the specified pattern
    regex = r'NODE_[0-9]+_length_[0-9]+_cov_[0-9]+\|\|full|NODE_[0-9]+_length_[0-9]+_cov_[0-9]+\|\|[0-9]+_partial'
    matches = re.findall(regex, text)
    return matches

def main(input_filename):
    # Debugging print to check the input file path
    print(f"Input filename: {input_filename}")
    
    # Read the input file
    with open(input_filename, 'r') as file:
        text = file.read()
    
    # Extract boundaries from the text
    boundaries = extract_boundaries(text)
    
    # Extract NODE details from the text
    nodes = extract_node_details(text)
    
    # Select the pair with the greatest range within the valid distance
    greatest_range = select_greatest_range(boundaries)
    
    # Create output filename in the same directory as the input file
    input_dir = os.path.dirname(input_filename)
    input_name = os.path.basename(input_filename)
    
    # Use the input filename to construct the output filename
    output_filename = os.path.join(input_dir, "{}_prophage_boundaries.tsv".format(os.path.splitext(input_name)[0]))
    
    # Write the result to the output file
    with open(output_filename, 'w') as file:
        if greatest_range and nodes:
            identifier, start, end = greatest_range
            file.write('{}\t{}\t{}\n'.format(nodes[0], start, end))
        else:
            file.write('No valid coordinate pairs found within the distance limit.\n')

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python extract_prophage_boundaries.py <input_filename>")
        sys.exit(1)
    
    input_filename = sys.argv[1]
    
    main(input_filename)
