from pathlib import Path

import numpy as np
import polars as pl
from joblib import load

loaded_model = load("trial2/output/phosphomannan_logreg.joblib")

embedding_df: pl.DataFrame = pl.read_parquet(
    "trial2/C_auris_embeddings_info.parquet",
)
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
phosphomannan_genes_df.write_csv(OUTPUT_DIR / "all_phosphomannan_genes.csv")
print(
    f"Percentage of phosphomannan genes: {len(phosphomannan_genes_df) / len(embedding_df) * 100:.2f}%"
)
