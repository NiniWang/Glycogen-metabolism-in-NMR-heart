from Bio import SeqIO
import pandas as pd
import re
import gzip

info_path = '{file_pathway}/longest_transcripts_info.csv'
raw_path = '{file_pathway}/pep.all.fa.gz'
final_path = '{file_pathway}/longest_protein.output'


file1_df = pd.read_csv(info_path)

index_elements = file1_df['original_index'].str.split().str[0]

with gzip.open(raw_path, 'rt') as raw_file:
    file2_records = SeqIO.to_dict(SeqIO.parse(raw_file, "fasta"))

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