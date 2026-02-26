#!/usr/bin/env bash

for file in *_1.fastq; do
  new_name="${file%_1.fastq}_R1.fastq"
  mv "$file" "$new_name"
done

for file in *_2.fastq; do
  new_name="${file%_2.fastq}_R2.fastq"
  mv "$file" "$new_name"
done