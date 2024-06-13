#step1 obtain the longest protein sequences for each of the two species separately.
python ExtractLongestProtein.py -info {file_pathway}/{species}_longest_transcripts_info.csv -raw {file_pathway}/{species}_pep.all.fa.gz -out {file_pathway}/{species}_longest_protein.fa

#step2 blastp
makeblastdb -in {species}_longest_protein.fa -dbtype prot -out {species}_longest_protein_db

blastp -query {file_pathway}/MM_longest_protein.fa -db {file_pathway}/NMR_longest_protein_db -out {file_pathway}/MMtoNMR_protein_blastp.output -evalue 1e-05  -max_target_seqs 1 -outfmt 6

blastp -query {file_pathway}/NMR_longest_protein.fa -db {file_pathway}/MM_longest_protein_db -out {file_pathway}/NMRtoMM_protein_blastp.output -evalue 1e-05  -max_target_seqs 1 -outfmt 6

#step3 merge BLASTP result files and find reciprocal hits
python find_reciprocal_hits.py MMtoNMR_protein_blastp.output NMRtoMM_protein_blastp.output MM_NMR_reciprocal_hits.output

#step4 Based on the results of reciprocal hits, merge the IDs of different species.
#Extract the ID information for each species.
python ExInfoFExInfoFromPep.py {file_pathway}/{species}_pep.all.fa.gz {species}_proteinID_geneID_transID.output

#Merge IDs based on protein ID correspondence
python mergeIDs.py {species1}_proteinID_geneID_transID.output {species2}_proteinID_geneID_transID.output {species1}_{species2}_reciprocal_hits.output {species1}_{species2}_allIDs_merge.output

