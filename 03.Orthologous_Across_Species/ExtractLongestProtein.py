import argparse
from Bio import SeqIO
import pandas as pd
import gzip

def process_files(info_path, raw_path, final_path):
    # Read the CSV file containing transcript information
    file1_df = pd.read_csv(info_path)

    # Extract the original index elements
    index_elements = file1_df['original_index'].str.split().str[0]

    # Parse the gzipped FASTA file
    with gzip.open(raw_path, 'rt') as raw_file:
        file2_records = SeqIO.to_dict(SeqIO.parse(raw_file, "fasta"))

    # Write the matched sequences to the output file
    with open(final_path, 'w') as output_file:
        for index_element in index_elements:
            for seq_id, seq_record in file2_records.items():
                protein_id = seq_id
                info = seq_record.description
                if index_element in info:
                    target_sequences = file2_records[protein_id]
                    SeqIO.write(target_sequences, output_file, "fasta")
                else:
                    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process files to extract specific protein sequences.")
    parser.add_argument("-info", help="Path to the CSV file containing transcript information", required=True)
    parser.add_argument("-raw", help="Path to the gzipped FASTA file", required=True)
    parser.add_argument("-out", help="Path to the output file", required=True)
    args = parser.parse_args()

    process_files(args.info, args.raw, args.out)