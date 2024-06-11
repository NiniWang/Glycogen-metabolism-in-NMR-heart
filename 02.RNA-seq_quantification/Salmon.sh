#step1
#build index
salmon index -t {species}_longest_trans_final.fa -i {species}_longest_trans_final_salmonindex -p 20

#step2
#quantification. Use 
sample=$1 #sample name
salmon quant -i {file_pathway}/{salmonindex} -l A -g {file_pathway}/{gtf_file} -1 {file_pathway}/${sample}_1_val_1.fq.gz -2 {file_pathway}/${sample}_2_val_2.fq.gz  -p 8 -o ${sample}

#step3
#merge results from same species
python merge_info.py