#!/usr/bin/env bash

#Functional Prediction 
#Run PICRUSt2 on ASV sequences and abundance table
picrust2_pipeline.py \
  -s 3-results/asv_sequences.fasta \
  -i 3-results/otu_table.biom \
  -o 3-results/picrust2_out_pipeline \
  -p 1

#Add KO descriptions
add_descriptions.py \
  -i 3-results/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
  -m KO \
  -o 3-results/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv

#Add MetaCyc pathway descriptions
add_descriptions.py \
  -i 3-results/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz \
  -m METACYC \
  -o 3-results/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv