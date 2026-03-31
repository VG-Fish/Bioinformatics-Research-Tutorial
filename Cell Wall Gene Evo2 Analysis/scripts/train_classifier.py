import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from joblib import dump
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import cross_val_score

train_embeddings_path: Path = Path("trial2/datasets/training_embeddings_info.parquet")
test_embeddings_path: Path = Path("trial2/datasets/test_embeddings_info.parquet")

if not train_embeddings_path.exists():
    raise FileNotFoundError(f"{train_embeddings_path} not found in working directory")
if not test_embeddings_path.exists():
    raise FileNotFoundError(f"{test_embeddings_path} not found in working directory")

OUTPUT_DIR: Path = Path("trial2/output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR / "figures", exist_ok=True)

print("Loading training data...")
train_df: pl.DataFrame = pl.read_parquet(train_embeddings_path)
X_train: np.ndarray = np.stack(train_df["Gene Embedding"].to_list())
y_train: np.ndarray = np.array(train_df["Is Phosphomannan Gene"].to_list(), dtype=int)
print(
    "Training data - X shape:",
    X_train.shape,
    "y shape:",
    y_train.shape,
    "positives:",
    y_train.sum(),
)

print("Loading test data...")
test_df: pl.DataFrame = pl.read_parquet(test_embeddings_path)
X_test: np.ndarray = np.stack(test_df["Gene Embedding"].to_list())
y_test: np.ndarray = np.array(test_df["Is Phosphomannan Gene"].to_list(), dtype=int)
print(
    "Test data - X shape:",
    X_test.shape,
    "y shape:",
    y_test.shape,
    "positives:",
    y_test.sum(),
)

print("\nTraining logistic regression model...")
logistic_regression_model: LogisticRegression = LogisticRegression(
    max_iter=1000,
    class_weight="balanced",
    random_state=42,
)
logistic_regression_model.fit(X_train, y_train)

output_model_path: Path = OUTPUT_DIR / "phosphomannan_logreg.joblib"
dump(logistic_regression_model, output_model_path)
print(f"Model saved to {output_model_path}")

print("\nEvaluating on test set (C. auris)...")
y_prob: np.ndarray = logistic_regression_model.predict_proba(X_test)[:, 1]
auroc: float = roc_auc_score(y_test, y_prob)  # pyright: ignore[reportAssignmentType]
print(f"Test AUROC: {auroc:.3f}")

print("Evaluating on training set (other species)...")
y_train_prob: np.ndarray = logistic_regression_model.predict_proba(X_train)[:, 1]
train_auroc: float = roc_auc_score(y_train, y_train_prob)  # pyright: ignore[reportAssignmentType]
print(f"Training AUROC: {train_auroc:.3f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

fpr_test, tpr_test, _ = roc_curve(y_test, y_prob)
ax1.plot(
    fpr_test, tpr_test, label=f"Test AUROC = {auroc:.3f}", color="red", linewidth=2
)
ax1.plot([0, 1], [0, 1], "k--", alpha=0.5)
ax1.set_xlabel("False Positive Rate")
ax1.set_ylabel("True Positive Rate")
ax1.set_title("Test Set ROC (C. auris)")
ax1.legend(loc="lower right")
ax1.grid(True, alpha=0.3)

fpr_train, tpr_train, _ = roc_curve(y_train, y_train_prob)
ax2.plot(
    fpr_train,
    tpr_train,
    label=f"Training AUROC = {train_auroc:.3f}",
    color="blue",
    linewidth=2,
)
ax2.plot([0, 1], [0, 1], "k--", alpha=0.5)
ax2.set_xlabel("False Positive Rate")
ax2.set_ylabel("True Positive Rate")
ax2.set_title("Training Set ROC (Other Species)")
ax2.legend(loc="lower right")
ax2.grid(True, alpha=0.3)

plt.tight_layout()
roc_path: Path = OUTPUT_DIR / "figures/roc_curves_comparison.png"
plt.savefig(roc_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"ROC curves saved to: {roc_path}")

thresholds: tuple[float, float, float, float] = (0.3, 0.5, 0.7, 0.8)
print("\nDetailed evaluation on test set:")
print("Threshold\tAccuracy\tPrecision\tRecall\tF1")
print("-" * 50)

best_threshold: float = 0.5
best_f1: float = 0

for threshold in thresholds:
    y_pred: np.ndarray = (y_prob >= threshold).astype(int)

    acc: float = accuracy_score(y_test, y_pred)

    tp: int = np.sum((y_pred == 1) & (y_test == 1))
    fp: int = np.sum((y_pred == 1) & (y_test == 0))
    fn: int = np.sum((y_pred == 0) & (y_test == 1))

    precision: float = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall: float = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1: float = (
        2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    )

    print(
        f"{threshold:.1f}\t\t{acc:.3f}\t\t{precision:.3f}\t\t{recall:.3f}\t\t{f1:.3f}"
    )

    if f1 >= best_f1:
        best_f1 = f1
        best_threshold = threshold

print(f"\nBest threshold based on F1-score: {best_threshold} (F1 = {best_f1:.3f})")

y_pred_final: np.ndarray = (y_prob >= best_threshold).astype(int)
acc_final: float = accuracy_score(y_test, y_pred_final)
cm: np.ndarray = confusion_matrix(y_test, y_pred_final)

print(f"\nFinal Results (threshold = {best_threshold}):")
print(f"Test Accuracy: {acc_final * 100:.2f}%")
print("Confusion Matrix:")
print(cm)
print(f"True Negatives: {cm[0, 0]}, False Positives: {cm[0, 1]}")
print(f"False Negatives: {cm[1, 0]}, True Positives: {cm[1, 1]}")

print("\nClassification Report:")
print(
    classification_report(
        y_test, y_pred_final, target_names=["Non-Phosphomannan", "Phosphomannan"]
    )
)

pred_df: pl.DataFrame = pl.DataFrame(
    {
        "Gene_Name": test_df["Gene Name"].to_list(),
        "y_true": y_test,
        "y_prob": y_prob,
        "y_pred": y_pred_final,
        "is_correct": (y_test == y_pred_final).astype(int),
    }
)

pred_path: Path = OUTPUT_DIR / "test_predictions.csv"
pred_df.write_csv(pred_path)
print(f"Test predictions saved to: {pred_path}")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print("Cross-species validation results:")
print(
    "Training species: C_albicans, C_dubliniensis, C_tropicalis, C_parapsilosis, C_lusitaniae, C_glabrata"
)
print("Test species: C_auris")
print(f"Test AUROC: {auroc:.3f}")
print(f"Test Accuracy: {acc_final:.3f}")
print(f"Best threshold: {best_threshold}")
print("Model ready for prediction on unlabeled C. auris genes")

cv_scores = cross_val_score(
    logistic_regression_model, X_train, y_train, cv=5, scoring="roc_auc"
)
print(f"Cross-validation AUROC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

embeddings_2d = TSNE(n_components=2).fit_transform(X_test)
plt.scatter(
    embeddings_2d[y_test == 0, 0],
    embeddings_2d[y_test == 0, 1],
    alpha=0.6,
    label="Negative",
)
plt.scatter(
    embeddings_2d[y_test == 1, 0],
    embeddings_2d[y_test == 1, 1],
    alpha=0.8,
    label="Positive",
)
plt.legend()
plt.title("Gene Embeddings Visualization")
plt.savefig(OUTPUT_DIR / "figures" / "gene_embeddings_visualization.png")
