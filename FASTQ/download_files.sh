#!/usr/bin/env bash

# This script requires GNU parallel (https://www.gnu.org/software/parallel/) to be installed
input_prefix="SRR"
declare -A experiment_groups=(
  [Experiment1]="$(printf '%s ' ${input_prefix}{11841493..11841495})"
  [Experiment2]="$(printf '%s ' ${input_prefix}{11841490..11841492})"
  [Experiment3]="$(printf '%s ' ${input_prefix}{11841487..11841489})"
  [Experiment4]="$(printf '%s ' ${input_prefix}{11841484..11841486})"
)

process_file() {
  read -r key val replicate_number <<< "$1"
  save_file="${key}_Replicate${replicate_number}.fastq"
  printf "Creating %s...\n" "$save_file"
  
  # Create temp directory and use trap to ensure that fasterq create directories are always cleaned up
  temp_job_dir=$(mktemp -d -t fasterq_temp_XXXXXX)
  trap 'rm -rf "${temp_job_dir}"' EXIT ERR INT

  # Quoting variables can prevent unexpected word-splitting errors
  fasterq-dump "$val" --split-files --threads 2 -t "$temp_job_dir" -o "$save_file" --mem 500MB

  purge
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