#!/usr/bin/env bash

#Download SRA runs, split F/R and convert to fastq, compress the output
while read -r srr; do
    [[ -z "$srr" ]] && continue
    fasterq-dump --split-files -O fastq "$srr" || {
        echo "Failed: $srr"
        continue
    }
    gzip fastq/"${srr}"_*.fastq
done < SRR_Acc_List.txt
