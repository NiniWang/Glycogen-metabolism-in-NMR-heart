
import argparse

def load_results_to_dict(file):
    results_dict = {}
    with open(file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            protein_id = parts[header.index('protein_id')]
            results_dict[protein_id] = line.strip()
    return results_dict, header

def read_mapping_file(mapping_file):
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                mapping[parts[0]] = parts[1]
    return mapping

def merge_files(species1_dict, species2_dict, mapping, output_file, header1, header2):
    with open(output_file, 'w') as f_out:
        # Write header
        f_out.write('\t'.join(header1) + '\t' + '\t'.join(header2) + '\n')
        
        for species1_protein_id, species2_protein_id in mapping.items():
            if species1_protein_id in species1_dict and species2_protein_id in species2_dict:
                species1_line = species1_dict[species1_protein_id]
                species2_line = species2_dict[species2_protein_id]
                f_out.write(species1_line + '\t' + species2_line + '\n')

def main():
    parser = argparse.ArgumentParser(description="Merge result files based on protein ID correspondence")
    parser.add_argument('species1_file', type=str, help="Input file for species 1")
    parser.add_argument('species2_file', type=str, help="Input file for species 2")
    parser.add_argument('mapping_file', type=str, help="Protein ID mapping file")
    parser.add_argument('output_file', type=str, help="Output file for merged results")
    
    args = parser.parse_args()
    
    # Load the results files into dictionaries
    species1_dict, header1 = load_results_to_dict(args.species1_file)
    species2_dict, header2 = load_results_to_dict(args.species2_file)
    
    # Read the protein ID mapping file
    mapping = read_mapping_file(args.mapping_file)
    
    # Merge the files and output the results
    merge_files(species1_dict, species2_dict, mapping, args.output_file, header1, header2)

if __name__ == "__main__":
    main()

