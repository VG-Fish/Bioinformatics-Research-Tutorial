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
export output_dir

# Make output_dir path absolute
output_dir=$(cd "$(dirname "$output_dir")" && pwd)/"$(basename "$output_dir")"
mkdir -p "$output_dir"
export output_dir


# Fix script to be more flexible in choosing an adapter
cd "$input_dir"
# sed is a stream editor, in this command, the format is s/old/new/, meaning substitute old with new.
ls *_1.fastq \
  | sed 's/_1\.fastq$//' \
  | parallel --progress -j 2 '
      BASE=$(basename {1})

      trimmomatic PE -threads 8 -phred33 \
        "{1}_1.fastq" "{1}_2.fastq" \
        "$output_dir/${BASE}_R1.paired.fastq" \
        "$output_dir/${BASE}_R1.unpaired.fastq" \
        "$output_dir/${BASE}_R2.paired.fastq" \
        "$output_dir/${BASE}_R2.unpaired.fastq" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
      
      rm -rf "{1}_1.fastq" "{1}_2.fastq"
  '