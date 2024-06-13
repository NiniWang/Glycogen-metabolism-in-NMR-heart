import re
import argparse
import gzip

def extract_ids_from_fasta(fasta_file, output_file):
    with gzip.open(fasta_file, 'rt') as f_in, open(output_file, 'w') as f_out:
        # write title line to output file
        f_out.write("protein_id\tgene_id\ttranscript_id\n")
        
        for line in f_in:
            if line.startswith('>'):
                # Use regex to extract IDs from the annotation line
                protein_id_match = re.search(r'^>(\S+)', line)
                gene_id_match = re.search(r'gene:(\S+)', line)
                transcript_id_match = re.search(r'transcript:(\S+)', line)
                
                if protein_id_match and gene_id_match and transcript_id_match:
                    protein_id = protein_id_match.group(1)
                    gene_id = gene_id_match.group(1)
                    transcript_id = transcript_id_match.group(1)
                    
                    # Write the extracted IDs to the output file
                    f_out.write(f"{protein_id}\t{gene_id}\t{transcript_id}\n")

def main():
    parser = argparse.ArgumentParser(description="Extract IDs from gzipped FASTA file")
    parser.add_argument('input_file', type=str, help="Input gzipped FASTA file")
    parser.add_argument('output_file', type=str, help="Output file to save extracted IDs")
    
    args = parser.parse_args()
    
    extract_ids_from_fasta(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

