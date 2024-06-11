#!/bin/bash
# extract the longest transcript from a genome annotation file.

python GetLongestTransFromGTF_NMR.py -g Heterocephalus_glaber_female.Naked_mole-rat_maternal.110.chr.gtf.gz -fa Heterocephalus_glaber_female.Naked_mole-rat_maternal.dna.toplevel.fa -o NMR_longest_trans.fa
