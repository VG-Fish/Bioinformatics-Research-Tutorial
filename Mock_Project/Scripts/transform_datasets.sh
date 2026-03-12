#!/usr/bin/env bash

genome_dir=""
feature_genome_file=""
input_dir=""
output_dir=""

print_usage() {
  echo "Usage: $0 -g genome_dir -f feature_genome_file -i input_dir -o output_dir" >&2
}

while getopts "g:f:i:o:" opt; do
  case "$opt" in
    g) genome_dir="$OPTARG" ;;
    f) feature_genome_file="$OPTARG" ;;
    i) input_dir="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    # The below case is a wildcard case that matches everything
    \?) print_usage; exit 1 ;;
  esac
done

# -z checks if the string is empty
if [[ -z "$genome_dir" || -z "$feature_genome_file" || -z "$input_dir" || -z "$output_dir" ]]; then
  # >&2 redirects output to stderr
  print_usage
  exit 1
fi

# Export these variables to ensure that parallel subshells can access them
export genome_dir
export feature_genome_file
export output_dir

mkdir -p "$output_dir"

run_hisat_and_samtools() {
  IFS="|" read -r R1 R2 sample_count type <<< "$1"

  output_sam_file="${output_dir}/${type}_Sample${sample_count}_aligned.sam"
  output_bam_file="${output_sam_file%.sam}.bam"
  output_sorted_bam_file="${output_bam_file%_aligned.bam}_sorted_aligned.bam"

  hisat2 \
    -x "$genome_dir" -p 3 \
    -1 "$R1" -2 "$R2" \
    -S "$output_sam_file" \

  samtools view -S -b "$output_sam_file" > "$output_bam_file"
  rm -rf "$output_sam_file"

  samtools sort "$output_bam_file" -o "$output_sorted_bam_file"
  rm -rf "$output_bam_file"

  samtools index "$output_sorted_bam_file"

  echo "\n\n\nRunning featureCounts: $output_sorted_bam_file...\n"

  featureCounts \
    -p --countReadPairs \
    -a "$feature_genome_file" \
    -o "$output_dir/${type}_Sample${sample_count}_counts.txt" \
    "$output_sorted_bam_file"
  # rm -rf "$output_sorted_bam_file"
}
export -f run_hisat_and_samtools

sample_count=1
for file in "$input_dir"/*_R1.paired.fastq; do
  R1="$file"
  R2="${file%_R1.paired.fastq}_R2.paired.fastq"

  test="ACOD"
  control="WTBM"
  type=""

  if [[ "$R1" == *"$test"* ]]; then
    type="$test"
  elif [[ "$R1" == *"$control"* ]]; then
    type="$control"
  else
    type="UNKNOWN"
  fi

  printf "%s|%s|%s|%s\n" "$R1" "$R2" "$sample_count" "$type"
  sample_count=$((sample_count+1))
done | parallel --compress --progress -j 6 run_hisat_and_samtools {}