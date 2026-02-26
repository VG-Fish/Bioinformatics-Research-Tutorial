#!/usr/bin/env bash

output_dir="../Trimmed FASTQ"
mkdir -p "$output_dir"

# sed is a stream editor, in this command, the format is s/old/new/, meaning substitute old with new.
ls *_R1.fastq \
  | sed 's/_R1\.fastq$//' \
  | parallel --progress -j 2 '
      R1={}_R1
      R2={}_R2
      trimmomatic PE -threads 16 -phred33 "${R1}.fastq" "${R2}.fastq" \
        "'"$output_dir"'/${R1}.paired.fastq" \
        "'"$output_dir"'/${R1}.unpaired.fastq" \
        "'"$output_dir"'/${R2}.paired.fastq" \
        "'"$output_dir"'/${R2}.unpaired.fastq" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  '