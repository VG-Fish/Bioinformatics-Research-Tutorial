#!/usr/bin/env bash

input_dir=""
output_dir=""

print_usage() {
  echo "Usage: $0 -i input_dir -o output_dir" >&2
}

while getopts "i:o:" opt; do
  case "$opt" in
    i) input_dir="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    # The below case is a wildcard case that matches everything
    \?) print_usage; exit 1 ;;
  esac
done

# -z checks if the string is empty
if [[ -z "$input_dir" || -z "$output_dir" ]]; then
  # >&2 redirects output to stderr
  print_usage
  exit 1
fi

# -p causes mkdir to not error and create parents if needed
mkdir -p "$output_dir"

for file in "$input_dir"/*.fastq; do
  printf "%s|%s\n" "$output_dir" "$file"
done | parallel --bar -j 4 '
  IFS="|" read -r output_dir input_file <<< {}
  fastqc -o "$output_dir" "$input_file"
'