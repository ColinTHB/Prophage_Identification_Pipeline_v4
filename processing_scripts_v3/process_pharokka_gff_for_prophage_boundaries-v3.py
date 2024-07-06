import sys
import os

def find_closest_cds(target_gene, cds_list):
    closest_cds = None
    min_distance = float('inf')
    
    for cds in cds_list:
        distance = abs(target_gene[2] - cds[2])
        if distance < min_distance:
            min_distance = distance
            closest_cds = cds
    
    return closest_cds, min_distance

def calculate_boundaries(gene, closest_gene):
    start = min(gene[2], closest_gene[2])
    end = max(gene[3], closest_gene[3])
    return start, end

def process_gff(input_file, output_dir):
    terminase_large_subunit = set()
    portal_protein = set()
    integrase = []
    endolysin = []
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            
            attributes = parts[8].split(';')
            gene_id = None
            product = None
            
            for attr in attributes:
                key_value = attr.split('=')
                if len(key_value) == 2:
                    if key_value[0] == 'ID':
                        gene_id = key_value[1]
                    elif key_value[0] == 'product':
                        product = key_value[1]
            
            gene_info = (parts[0], gene_id, int(parts[3]), int(parts[4]))
            if product == 'terminase large subunit':
                terminase_large_subunit.add(gene_info)
            elif product == 'portal protein':
                portal_protein.add(gene_info)
            elif product == 'integrase':
                integrase.append(gene_info)
            elif product == 'endolysin':
                endolysin.append(gene_info)

    output_file = os.path.join(output_dir, os.path.basename(input_file).rsplit('.', 1)[0] + '_output.txt')
    with open(output_file, 'w') as out:
        terminase_count = 1
        for target_gene in terminase_large_subunit:
            closest_integrase, dist_integrase = find_closest_cds(target_gene, integrase)
            closest_endolysin, dist_endolysin = find_closest_cds(target_gene, endolysin)
            
            if closest_integrase:
                start_integrase, end_integrase = closest_integrase[2], closest_integrase[3]
            else:
                start_integrase, end_integrase = target_gene[2], target_gene[3]
            
            if closest_endolysin:
                start_endolysin, end_endolysin = closest_endolysin[2], closest_endolysin[3]
            else:
                start_endolysin, end_endolysin = target_gene[2], target_gene[3]
            
            out.write(f"Closest integrase for terminase {terminase_count}: {closest_integrase[0]} {closest_integrase[1]} {closest_integrase[2]} {closest_integrase[3]}\n")
            out.write(f"Closest endolysin for terminase {terminase_count}: {closest_endolysin[0]} {closest_endolysin[1]} {closest_endolysin[2]} {closest_endolysin[3]}\n")
            out.write(f"Boundaries for terminase {terminase_count}: {min(start_integrase, start_endolysin)} to {max(end_integrase, end_endolysin)}\n")
            
            terminase_count += 1
        
        portal_count = 1
        for target_gene in portal_protein:
            closest_integrase, dist_integrase = find_closest_cds(target_gene, integrase)
            closest_endolysin, dist_endolysin = find_closest_cds(target_gene, endolysin)
            
            if closest_integrase:
                start_integrase, end_integrase = closest_integrase[2], closest_integrase[3]
            else:
                start_integrase, end_integrase = target_gene[2], target_gene[3]
            
            if closest_endolysin:
                start_endolysin, end_endolysin = closest_endolysin[2], closest_endolysin[3]
            else:
                start_endolysin, end_endolysin = target_gene[2], target_gene[3]
            
            out.write(f"Closest integrase for portal protein {portal_count}: {closest_integrase[0]} {closest_integrase[1]} {closest_integrase[2]} {closest_integrase[3]}\n")
            out.write(f"Closest endolysin for portal protein {portal_count}: {closest_endolysin[0]} {closest_endolysin[1]} {closest_endolysin[2]} {closest_endolysin[3]}\n")
            out.write(f"Boundaries for portal protein {portal_count}: {min(start_integrase, start_endolysin)} to {max(end_integrase, end_endolysin)}\n")
            
            portal_count += 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_gff_file output_dir")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    process_gff(input_file, output_dir)
    print(f"Output written to {os.path.join(output_dir, os.path.basename(input_file).rsplit('.', 1)[0])}_output.txt")
