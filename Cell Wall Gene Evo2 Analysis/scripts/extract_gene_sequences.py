import random
from os import makedirs
from pathlib import Path

import polars as pl
from Bio import SeqIO

CDS_DATA_DIRECTORY: Path = Path("CDS Data")
OUTPUT_DIRECTORY: Path = Path("trial2/sequences_output")
if not OUTPUT_DIRECTORY.exists():
    makedirs(OUTPUT_DIRECTORY)

organisms: list[str] = [
    "C_albicans",
    "C_auris",
    "C_dubliniensis",
    "C_tropicalis",
    "C_parapsilosis",
    "C_lusitaniae",
    "C_glabrata",
]

PHOSPHOMANNAN_GENE_NAMES: str = "Phosphomannan Genes in Candida species.csv"
phosphomannan_genes_df: pl.DataFrame = pl.read_csv(PHOSPHOMANNAN_GENE_NAMES).select(
    organisms
)

random_gene_multiplier: float = 1.5

gene_sequences: list[dict[str, str | bool]] = []
for organism in phosphomannan_genes_df.columns:
    genes: list[str] = phosphomannan_genes_df[organism].drop_nulls().to_list()
    full_path: Path = next(CDS_DATA_DIRECTORY.glob(f"{organism}*.fasta"))

    id_to_seq: dict = {
        record.id: str(record.seq) for record in SeqIO.parse(full_path, "fasta")
    }

    num_genes: int = 0
    for gene in genes:
        sequence = id_to_seq.get(gene)
        if sequence:
            num_genes += 1
            gene_sequences.append(
                {"Gene Id": gene, "Sequence": sequence, "Phosphomannan Gene": True}
            )

    num_genes = int(num_genes * random_gene_multiplier)
    all_genes: list[str] = list(id_to_seq.keys())
    random_genes: list[str] = random.sample(all_genes, num_genes)
    for gene in random_genes:
        sequence = id_to_seq[gene]
        gene_sequences.append(
            {"Gene Id": gene, "Sequence": sequence, "Phosphomannan Gene": False}
        )

gene_sequences_df: pl.DataFrame = pl.DataFrame(
    gene_sequences,
    schema={
        "Name": pl.String,
        "Sequence": pl.String,
        "Phosphomannan Gene": pl.Boolean,
    },
)

gene_sequences_path: Path = OUTPUT_DIRECTORY / "gene_sequences.csv"
gene_sequences_df.write_csv(gene_sequences_path)
