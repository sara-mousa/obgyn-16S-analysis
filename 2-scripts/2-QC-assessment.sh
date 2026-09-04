#!/usr/bin/env bash

#Create directories for QC results
mkdir -p 3-results/1-qc/fastqc
mkdir -p 3-results/1-qc/multiqc

#Run FastQC on compressed FASTQ files
fastqc 1-data/fastq/*.fastq.gz -o 3-results/1-qc/fastqc

#Summarize QC results using MultiQC
multiqc 3-results/qc/fastqc -o 3-results/1-qc/multiqc
