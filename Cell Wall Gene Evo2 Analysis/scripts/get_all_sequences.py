from os import makedirs
from pathlib import Path

import polars as pl
from Bio import SeqIO

CDS_DATA_DIRECTORY: Path = Path("CDS Data")
full_path: Path = (
    CDS_DATA_DIRECTORY / "C_auris_B8441_version_s02-m01-r12_orf_coding.fasta"
)
OUTPUT_DIRECTORY: Path = Path("trial2")
makedirs(OUTPUT_DIRECTORY, exist_ok=True)

records = list(SeqIO.parse(full_path, "fasta"))
gene_sequences_df = pl.DataFrame(
    {
        "Name": [r.id for r in records],
        "Sequence": [str(r.seq) for r in records],
    }
)
gene_sequences_df.write_csv(OUTPUT_DIRECTORY / "C_auris_gene_sequences.csv")
