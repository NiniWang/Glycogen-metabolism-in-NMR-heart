# PRECALCULATED BLAST RESULTS
# COMMENT OUT THESE TWO LINES IF YOU WANT TO USE THE RESULTS FROM THE CELL ABOVE
fwd_out = open('MMtoNMR_protein_blastp.output')
rev_out = open('NMRtoMM_protein_blastp.output')

# Load the BLAST results into Pandas dataframes
fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
rev_results = pd.read_csv(rev_out, sep="\t", header=None)

# Add headers to forward and reverse results dataframes
headers = ["query", "subject", "identity", "coverage",
           "qlength", "slength", "alength",
           "bitscore", "E-value"]
fwd_results.columns = headers
rev_results.columns = headers

# Create a new column in both dataframes: normalised bitscore
fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

# Create query and subject coverage columns in both dataframes
fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
rev_results['qcov'] = rev_results.alength/rev_results.qlength
fwd_results['scov'] = fwd_results.alength/fwd_results.slength
rev_results['scov'] = rev_results.alength/rev_results.slength

# Clip maximum coverage values at 1.0
fwd_results['qcov'] = fwd_results['qcov'].clip_upper(1)
rev_results['qcov'] = rev_results['qcov'].clip_upper(1)
fwd_results['scov'] = fwd_results['scov'].clip_upper(1)
rev_results['scov'] = rev_results['scov'].clip_upper(1)

# Merge forward and reverse results
rbbh = pd.merge(fwd_results, rev_results[['query', 'subject']],
                left_on='subject', right_on='query',
                how='outer')

# Discard rows that are not RBH
rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]

# Group duplicate RBH rows, taking the maximum value in each column
rbbh = rbbh.groupby(['query_x', 'subject_x']).max()
