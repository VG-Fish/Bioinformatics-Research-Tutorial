from pathlib import Path

import polars as pl
from Bio import SeqIO

CDS_DATA_DIRECTORY: Path = Path("CDS Data")
full_path: Path = CDS_DATA_DIRECTORY / "C_auris_B11221_orf_coding.fasta"

records = list(SeqIO.parse(full_path, "fasta"))
gene_sequences_df = pl.DataFrame(
    {
        "Name": [r.id for r in records],
        "Sequence": [str(r.seq) for r in records],
    }
)
gene_sequences_df.write_csv("C_auris_gene_sequences.csv")
