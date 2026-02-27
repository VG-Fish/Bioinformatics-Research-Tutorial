#!/usr/bin/env bash

# This script requires GNU parallel (https://www.gnu.org/software/parallel/) to be installed
# Got SRR ids from: https://www.ncbi.nlm.nih.gov/sra/?term=SRP263458

output_dir=""

print_usage() {
  echo "Usage: $0 -o output_dir" >&2
}

while getopts "o:" opt; do
  case "$opt" in
    o) output_dir="$OPTARG" ;;
    # The below case is a wildcard case that matches everything
    \?) print_usage; exit 1 ;;
  esac
done

# -z checks if the string is empty
if [[ -z "$output_dir" ]]; then
  # >&2 redirects output to stderr
  print_usage
  exit 1
fi

# -p causes mkdir to not error and create parents if needed
export output_dir
mkdir -p "$output_dir"


input_prefix="SRR"
declare -A experiment_groups=(
  [Experiment1]="$(printf '%s ' ${input_prefix}{11841493..11841494})"
  [Experiment2]="$(printf '%s ' ${input_prefix}{11841490..11841491})"
  # [Experiment3]="$(printf '%s ' ${input_prefix}{11841487..11841489})"
  # [Experiment4]="$(printf '%s ' ${input_prefix}{11841484..11841486})"
)

process_file() {
  read -r key val replicate_number <<< "$1"
  save_file="$output_dir/${key}_Replicate${replicate_number}.fastq"
  printf "Creating %s...\n" "$save_file"
  
  # Create temp directory and use trap to ensure that fasterq create directories are always cleaned up
  temp_job_dir=$(mktemp -d -t fasterq_temp_XXXXXX)
  trap 'rm -rf "${temp_job_dir}"' EXIT ERR INT

  printf "1. Prefetching %s...\n" "$val"
  prefetch "$val" -O "$temp_job_dir"

  # Quoting variables can prevent unexpected word-splitting errors
  printf "2. Extracting %s to %s...\n" "$val" "$save_file"
  fasterq-dump "${temp_job_dir}/${val}/${val}.sra" \
    --split-files -x \
    --threads 2 --mem 500MB \
    -t "$temp_job_dir" -o "$save_file"
}
export -f process_file

for key in "${!experiment_groups[@]}"; do
  vals=${experiment_groups[$key]}
  replicate_number=1
  for v in $vals; do
    printf "%s %s %s\n" "$key" "$v" "$replicate_number"
    replicate_number=$((replicate_number+1))
  done
done | parallel --progress -j 4 process_file {}

# Renaming the files to have R1 & R2 at the end to be more clear
for file in "$output_dir"/*_1.fastq; do
  new_name="${file%_1.fastq}_R1.fastq"
  mv "$file" "$new_name"
done

for file in "$output_dir"/*_2.fastq; do
  new_name="${file%_2.fastq}_R2.fastq"
  mv "$file" "$new_name"
done