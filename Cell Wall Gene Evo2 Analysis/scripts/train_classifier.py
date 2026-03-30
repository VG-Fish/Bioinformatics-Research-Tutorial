import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from joblib import dump
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split

embeddings_path: Path = Path("trial2/training_embeddings_info.parquet")
if not embeddings_path.exists():
    raise FileNotFoundError(f"{embeddings_path} not found in working directory")

OUTPUT_DIR: Path = Path("trial2/output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR / "figures", exist_ok=True)

df: pl.DataFrame = pl.read_parquet(embeddings_path)
X: np.ndarray = np.stack(df["Gene Embedding"].to_list())
y: np.ndarray = np.array(df["Is Phophomannam Gene"].to_list(), dtype=int)
print("X shape:", X.shape, "y shape:", y.shape, "positives:", y.sum())

RANDOM_SEED: int = 42
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.3,
    random_state=RANDOM_SEED,
    stratify=y,
)

logistic_regression_model: LogisticRegression = LogisticRegression(
    max_iter=1000,
    class_weight="balanced",
)
logistic_regression_model.fit(X_train, y_train)
output_model_path: Path = OUTPUT_DIR / "phosphomannan_logreg.joblib"
dump(logistic_regression_model, output_model_path)
print(f"Model saved to {output_model_path}")

y_prob: np.ndarray = logistic_regression_model.predict_proba(X_test)[:, 1]
auroc: float = roc_auc_score(y_test, y_prob)  # pyright: ignore[reportAssignmentType]
print(f"Test AUROC: {auroc:.3f}")

fpr, tpr, thresholds = roc_curve(y_test, y_prob)

plt.figure(figsize=(4, 4))
plt.plot(fpr, tpr, label=f"AUROC = {auroc:.3f}")
plt.plot([0, 1], [0, 1], "k--", alpha=0.5)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Phosphomannan Gene Classifier ROC")
plt.legend(loc="lower right")
plt.tight_layout()

os.makedirs("figures", exist_ok=True)
roc_path: Path = OUTPUT_DIR / "figures/roc_curve.png"
plt.savefig(roc_path, dpi=150)
plt.close()

print("ROC curve saved to:", roc_path)

threshold: float = 0.5
y_pred: np.ndarray = (y_prob >= threshold).astype(int)

acc: float = accuracy_score(y_test, y_pred)
cm: np.ndarray = confusion_matrix(y_test, y_pred)
print(f"Accuracy: {acc * 100:.2f}%")
print("Confusion matrix:\n", cm)

pred_df: pl.DataFrame = pl.DataFrame(
    {"y_true": y_test, "y_prob": y_prob, "y_pred": y_pred}
)
os.makedirs("outputs", exist_ok=True)
pred_path: Path = OUTPUT_DIR / "test_predictions.csv"
pred_df.write_csv(pred_path)
print("Test predictions saved to:", pred_path)
