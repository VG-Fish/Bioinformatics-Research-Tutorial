from pathlib import Path

import numpy as np
import polars as pl
from joblib import load

loaded_model = load("trial2/output/phosphomannan_logreg.joblib")

PHOSPHOMANNAN_GENE_NAMES: str = "Phosphomannan Genes in Candida species.csv"
c_auris_phosphomannan_genes: list[str] = (
    pl.read_csv(PHOSPHOMANNAN_GENE_NAMES).select(pl.col("C_auris")).drop_nulls()
)["C_auris"].to_list()

embedding_df: pl.DataFrame = pl.read_parquet(
    "trial2/C_auris_embeddings_info.parquet",
).filter(pl.col("Gene Name").is_in(c_auris_phosphomannan_genes))
OUTPUT_DIR: Path = Path("trial2/output")

phosphomannan_genes_info: list[dict[str, str | float]] = []
# phosphomannan_gene_threshold: float = 0.5

for embedding_info in embedding_df.rows(named=True):
    X_new: np.ndarray = np.array(embedding_info["Gene Embedding"]).reshape(1, -1)
    y_prob: float = loaded_model.predict_proba(X_new)[0, 1]
    phosphomannan_genes_info.append(
        {"Gene Name": embedding_info["Gene Name"], "Model Confidence": y_prob}
    )

phosphomannan_genes_df: pl.DataFrame = pl.DataFrame(phosphomannan_genes_info).sort(
    by=pl.col("Model Confidence"), descending=True
)
phosphomannan_genes_df.write_csv(OUTPUT_DIR / "confirmed_phosphomannan_genes.csv")
print(
    f"Percentage of phosphomannan genes: {len(phosphomannan_genes_df) / len(embedding_df) * 100:.2f}%"
)
