#!/usr/bin/env bash

genome_dir=""
input_dir=""
output_dir=""

print_usage() {
  echo "Usage: $0 -g genome_dir -i input_dir -o output_dir" >&2
}

while getopts "g:i:o:" opt; do
  case "$opt" in
    g) genome_dir="$OPTARG" ;;
    i) input_dir="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    # The below case is a wildcard case that matches everything
    \?) print_usage; exit 1 ;;
  esac
done

# -z checks if the string is empty
if [[ -z "$genome_dir" || -z "$input_dir" || -z "$output_dir" ]]; then
  # >&2 redirects output to stderr
  print_usage
  exit 1
fi

# Export these variables to ensure that parallel subshells can access them
export genome_dir
export output_dir

mkdir -p "$output_dir"

run_hisat() {
  IFS="|" read -r R1 R2 sample_count <<< "$1"
  hisat2 \
    -x "$genome_dir" -p 3 \
    -1 "$R1" -2 "$R2" \
    -S "${output_dir}/Sample${sample_count}_aligned.sam"
}
export -f run_hisat

sample_count=1
for file in "$input_dir"/*_R1.paired.fastq; do
  R1="$file"
  R2="${file%_R1.paired.fastq}_R2.paired.fastq"
  printf "%s|%s|%s\n" "$R1" "$R2" "$sample_count"
  sample_count=$((sample_count+1))
done | parallel --progress -j 4 run_hisat {}