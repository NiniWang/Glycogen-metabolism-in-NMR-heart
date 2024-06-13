#step1 clean reads
sample=$1
trim_galore -q 20 --paired {file_pathway}/${sample}_1.fastq.gz {file_pathway}/${sample}_2.fastq.gz -o {file_pathway}/${sample}

#step2 build index
salmon index -t {species}_longest_trans.fa -i {species}_longest_trans_salmonindex -p 20

#step3 quantification 
sample=$1 #sample name
salmon quant -i {file_pathway}/{salmonindex} -l A -g {file_pathway}/{gtf_file} -1 {file_pathway}/${sample}/${sample}_1_val_1.fq.gz -2 {file_pathway}/${sample}/${sample}_2_val_2.fq.gz  -p 8 -o {file_pathway}/${sample}

#step4 merge results from same species
python merge_salmon_result.py -folders {sample1} {sample2} {sample3}  -col NumReads -out merged_{species}_NumReads.csv
python merge_salmon_result.py -folders {sample1} {sample2} {sample3}  -col TPM -out merged_{species}_TPM.csv
python merge_salmon_result.py -folders {sample1} {sample2} {sample3}  -col EffectiveLength -out merged_{species}_EffectiveLength.csv

#step5 merge NumReads/EffectiveLength/TPM results from different species into same matrix based on the result in 03.Orthologous_Across_Species
python MergeResultAcrossSpecies.py merged_{species1}_{NumReads/EffectiveLength/TPM}.csv merged_{species2}_{NumReads/EffectiveLength/TPM}.csv {species1}_{species2}_allIDs_merge.output merged_{species1}_{species2}_{NumReads/EffectiveLength/TPM}.output

#step6 identify DEGs
# see Identify_DEGs.R script