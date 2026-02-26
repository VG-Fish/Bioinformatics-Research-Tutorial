#!/usr/bin/env bash

output_dir=""

while getopts "o:" opt; do
  case "$opt" in
    o) output_dir="$OPTARG" ;;
    # The below case is a wildcard cases that matches everything
    \?) echo "Usage: $0 -o output_dir" >&2; exit 1 ;;
  esac
done

# -z checks if the string is empty
if [[ -z "$output_dir" ]]; then
  # >&2 redirects output to stderr
  echo "Usage: $0 -o output_dir" >&2
  exit 1
fi

# -p causes mkdir to not error and create parents if needed
mkdir -p "$output_dir"

# ::: is shell globbing
parallel --progress -j 4 fastqc -o "$output_dir" ::: FASTQ/*.fastq