import random
from os import makedirs
from pathlib import Path

import polars as pl
from Bio import SeqIO

CDS_DATA_DIRECTORY: Path = Path("CDS Data")
OUTPUT_DIRECTORY: Path = Path("trial2/sequences_output")
if not OUTPUT_DIRECTORY.exists():
    makedirs(OUTPUT_DIRECTORY)

train_organisms: list[str] = [
    "C_albicans",
    "C_dubliniensis",
    "C_tropicalis",
    "C_parapsilosis",
    "C_lusitaniae",
    "C_glabrata",
    "S_cerevisiae",
]
test_organism: str = "C_auris"

PHOSPHOMANNAN_GENE_NAMES: str = "Phosphomannan Genes in Candida species.csv"
phosphomannan_genes_df: pl.DataFrame = pl.read_csv(PHOSPHOMANNAN_GENE_NAMES)
random_gene_multiplier: float | int = 5


def get_gene_sequences_df(*organisms) -> pl.DataFrame:
    gene_sequences: list[dict[str, str | bool]] = []

    for organism in organisms:
        genes: list[str] = phosphomannan_genes_df[organism].drop_nulls().to_list()
        full_path: Path = next(CDS_DATA_DIRECTORY.glob(f"{organism}*.fasta"))

        id_to_seq: dict = {
            record.id: str(record.seq) for record in SeqIO.parse(full_path, "fasta")
        }

        num_phosphomannan_genes: int = 0
        for gene in genes:
            sequence = id_to_seq.get(gene)
            if sequence:
                num_phosphomannan_genes += 1
                gene_sequences.append(
                    {
                        "Name": gene,
                        "Sequence": sequence,
                        "Is Phosphomannan Gene": True,
                    }
                )

        num_random_genes = int(num_phosphomannan_genes * random_gene_multiplier)
        all_genes: list[str] = list(id_to_seq.keys())

        num_random_genes = min(num_random_genes, len(all_genes))
        random_genes: list[str] = random.sample(all_genes, num_random_genes)

        for gene in random_genes:
            sequence = id_to_seq[gene]
            gene_sequences.append(
                {
                    "Name": gene,
                    "Sequence": sequence,
                    "Is Phosphomannan Gene": False,
                }
            )

        print(
            f"{organism}: {num_phosphomannan_genes} phosphomannan genes, {num_random_genes} random genes"
        )

    return pl.DataFrame(
        gene_sequences,
        schema={
            "Name": pl.String,
            "Sequence": pl.String,
            "Is Phosphomannan Gene": pl.Boolean,
        },
    )


print("Generating training data...")
train_gene_sequences_df: pl.DataFrame = get_gene_sequences_df(*train_organisms)
train_gene_sequences_path: Path = OUTPUT_DIRECTORY / "train_gene_sequences.csv"
train_gene_sequences_df.write_csv(train_gene_sequences_path)
print(f"Training data saved to: {train_gene_sequences_path}")

print("Generating test data...")
test_gene_sequences_df: pl.DataFrame = get_gene_sequences_df(test_organism)
test_gene_sequences_path: Path = OUTPUT_DIRECTORY / "test_gene_sequences.csv"
test_gene_sequences_df.write_csv(test_gene_sequences_path)
print(f"Test data saved to: {test_gene_sequences_path}")

print("Summary:")
print(f"Training set: {len(train_gene_sequences_df)} sequences")
print(f"Test set: {len(test_gene_sequences_df)} sequences")
print(
    f"Training phosphomannan genes: {train_gene_sequences_df['Is Phosphomannan Gene'].sum()}"
)
print(
    f"Test phosphomannan genes: {test_gene_sequences_df['Is Phosphomannan Gene'].sum()}"
)
