import argparse

def parse_blast_output(file):
    """
    Parse BLAST output file in format 6 and return a dictionary with query_id as key and subject_id as value.
    """
    blast_dict = {}
    with open(file, 'r') as f:
        for line in f:
            columns = line.strip().split()
            query_id = columns[0]
            subject_id = columns[1]
            blast_dict[query_id] = subject_id
    return blast_dict

def find_reciprocal_hits(file1, file2):
    """
    Find reciprocal hits between two BLAST output files.
    """
    blast1 = parse_blast_output(file1)
    blast2 = parse_blast_output(file2)
    
    reciprocal_hits = []
    
    for query_id in blast1:
        subject_id = blast1[query_id]
        if subject_id in blast2 and blast2[subject_id] == query_id:
            reciprocal_hits.append((query_id, subject_id))
    
    return reciprocal_hits

def main():
    parser = argparse.ArgumentParser(description="Find reciprocal hits between two BLAST output files.")
    parser.add_argument('file1', type=str, help="First BLAST output file.")
    parser.add_argument('file2', type=str, help="Second BLAST output file.")
    parser.add_argument('output', type=str, help="Output file for reciprocal hits.")
    
    args = parser.parse_args()
    
    file1 = args.file1
    file2 = args.file2
    output = args.output
    
    reciprocal_hits = find_reciprocal_hits(file1, file2)
    
    with open(output, 'w') as f:
        for query_id, subject_id in reciprocal_hits:
            f.write(f"{query_id}\t{subject_id}\n")

if __name__ == "__main__":
    main()
