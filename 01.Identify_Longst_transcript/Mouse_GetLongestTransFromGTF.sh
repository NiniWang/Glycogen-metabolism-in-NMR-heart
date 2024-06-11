#!/bin/bash
# extract the longest transcript from a genome annotation file.

python GetLongestTransFromGTF_Mouse.py -g Mus_musculus.GRCm39.110.chr.gtf.gz -fa Mus_musculus.GRCm39.dna.toplevel.fa -o MM_longest_trans.fa
