# Corrected re-run of scripts/03_biomarker_shap.py, on the padj-filtered
# feature matrix. Same RF / Keras NN / SHAP logic as the original, plus:
#   (a) tf.random.set_seed(42) for NN reproducibility
#   (b) classification_report + accuracy_score for RF -> rf_metrics.txt
#   (c) model.evaluate() for NN -> nn_metrics.txt
#   (d) Keras training history saved -> nn_training_history.csv
#   (e) same RF importance ranking + SHAP top-20/20 selection logic

import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shap

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
from sklearn.preprocessing import StandardScaler

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout

# CHANGE (a): seed the Keras/TensorFlow model, in addition to the existing
# sklearn seeds (train_test_split random_state=42, RandomForest random_state=52)
tf.random.set_seed(42)

PROJECT_DIR = "/Users/nishantthalwal/Documents/Data ALL/TCGA_LUAD"
DATA_PATH = os.path.join(PROJECT_DIR, "publication_v1", "data", "FEATURE_DATA_MATRIX.csv")
RESULTS_DIR = os.path.join(PROJECT_DIR, "publication_v1", "results")
FIG_DIR = os.path.join(PROJECT_DIR, "publication_v1", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# Load the data
data = pd.read_csv(DATA_PATH, index_col=0)
print("Data shape:", data.shape)

X = data.drop('Label', axis=1)
y = data['Label']

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
X_scaled_df = pd.DataFrame(X_scaled, columns=X.columns)

X_train, X_test, y_train, y_test = train_test_split(X_scaled_df, y, test_size=0.2, random_state=42)

# ==========================================
# 1. RANDOM FOREST MODEL & FEATURE IMPORTANCE
# ==========================================
rf_model = RandomForestClassifier(n_estimators=100, random_state=52)
rf_model.fit(X_train, y_train)

y_pred = rf_model.predict(X_test)

# CHANGE (b): RF performance metrics, saved (never computed in the original)
rf_report = classification_report(y_test, y_pred)
rf_accuracy = accuracy_score(y_test, y_pred)
with open(os.path.join(RESULTS_DIR, "rf_metrics.txt"), "w") as f:
    f.write("Random Forest test-set performance\n")
    f.write("===================================\n")
    f.write(f"Accuracy: {rf_accuracy:.6f}\n\n")
    f.write("Classification report:\n")
    f.write(rf_report)
print("RF accuracy:", rf_accuracy)

feature_importances = rf_model.feature_importances_
feature_names = X.columns

importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': feature_importances})
importance_df = importance_df.sort_values(by='Importance', ascending=False)

plt.figure(figsize=(12, 8))
plt.barh(importance_df['Feature'][:20], importance_df['Importance'][:20])
plt.xlabel('Importance')
plt.ylabel('Feature')
plt.title('Top 20 Genes : RandomForest')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.gca().invert_yaxis()
plt.savefig(os.path.join(FIG_DIR, "Top_MARKERS_RandomForest.png"), bbox_inches="tight")
plt.close()

# ==========================================
# 2. DEEP LEARNING MODEL & SHAP EXPLANATION
# ==========================================
model = Sequential([
    Dense(128, input_dim=X_train.shape[1], activation='relu'),
    Dropout(0.6),
    Dense(64, activation='relu'),
    Dropout(0.6),
    Dense(32, activation='relu'),
    Dropout(0.6),
    Dense(16, activation='relu'),
    Dropout(0.6),
    Dense(1, activation='sigmoid')
])
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
history = model.fit(X_train, y_train, epochs=50, batch_size=32, validation_split=0.2, verbose=1)

# CHANGE (d): save training history (never saved in the original)
history_df = pd.DataFrame(history.history)
history_df.insert(0, "epoch", range(1, len(history_df) + 1))
history_df.to_csv(os.path.join(RESULTS_DIR, "nn_training_history.csv"), index=False)

# CHANGE (c): evaluate on the held-out test set (never called in the original)
test_loss, test_accuracy = model.evaluate(X_test, y_test, verbose=0)
with open(os.path.join(RESULTS_DIR, "nn_metrics.txt"), "w") as f:
    f.write("Keras Neural Network test-set performance\n")
    f.write("==========================================\n")
    f.write(f"Test loss: {test_loss:.6f}\n")
    f.write(f"Test accuracy: {test_accuracy:.6f}\n")
print("NN test loss:", test_loss, " test accuracy:", test_accuracy)

X_train_np = X_train.values if isinstance(X_train, pd.DataFrame) else np.array(X_train)
X_test_np = X_test.values if isinstance(X_test, pd.DataFrame) else np.array(X_test)

explainer = shap.Explainer(model, X_train_np)
shap_values = explainer(X_test_np, max_evals=1000)

if hasattr(shap_values, "values"):
    shap_values_matrix = shap_values.values
else:
    shap_values_matrix = np.array(shap_values)

if len(shap_values_matrix.shape) == 3:
    shap_values_matrix = np.squeeze(shap_values_matrix, axis=-1)
elif shap_values_matrix.shape[0] != X_test_np.shape[0] and shap_values_matrix.shape[1] == X_test_np.shape[0]:
    shap_values_matrix = shap_values_matrix.T

shap_values_class = shap_values_matrix

y_test_aligned = y_test.reset_index(drop=True)
X_test_aligned = X_test.reset_index(drop=True)

mask_class_0 = (y_test_aligned == 0)
mask_class_1 = (y_test_aligned == 1)

shap_values_class_1 = shap_values_class[mask_class_1.values]
mean_signed_class_1 = shap_values_class_1.mean(0)  # signed, not abs: sign = direction of push
SHAP_CLASS_1 = pd.DataFrame({
    'Features': X_test_aligned.columns.tolist(),
    'Mean Signed SHAP': mean_signed_class_1
})

shap_values_class_0 = shap_values_class[mask_class_0.values]
mean_signed_class_0 = shap_values_class_0.mean(0)
SHAP_CLASS_0 = pd.DataFrame({
    'Features': X_test_aligned.columns.tolist(),
    'Mean Signed SHAP': mean_signed_class_0
})

# ==========================================
# 3. BIOMARKER SELECTION LOGIC
# ==========================================
# CHANGE: rank by signed mean SHAP contribution per class, instead of the
# original's abs-magnitude ranking + cross-class exclusivity swap. The
# original approach broke down on a smaller (110-gene) candidate pool: the
# tumor- and normal-subset importance rankings are highly correlated
# (Pearson r=0.87 here) since both come from one binary classifier's
# explanation, so the add-then-evict exclusivity loop emptied both lists
# before reaching 20 markers each (see investigation in this session).
#
# This replacement asks the more direct question: among tumor-labeled test
# samples, does this gene's SHAP value push the prediction toward tumor
# (positive) or away from it (negative)? TUMOR_MARKER = genes with the
# strongest positive mean signed SHAP among tumor samples. NORMAL_MARKER =
# genes with the strongest negative mean signed SHAP among normal samples.
# A gene can appear in both lists (no forced exclusivity) -- that is a
# meaningful outcome (a gene consistently discriminative in both
# directions), not an error.
TUMOR_MARKER = (
    SHAP_CLASS_1[SHAP_CLASS_1['Mean Signed SHAP'] > 0]
    .sort_values('Mean Signed SHAP', ascending=False)['Features']
    .head(20)
    .tolist()
)
NORMAL_MARKER = (
    SHAP_CLASS_0[SHAP_CLASS_0['Mean Signed SHAP'] < 0]
    .sort_values('Mean Signed SHAP', ascending=True)['Features']
    .head(20)
    .tolist()
)
print(f"TUMOR_MARKER selected: {len(TUMOR_MARKER)} (of up to 20)")
print(f"NORMAL_MARKER selected: {len(NORMAL_MARKER)} (of up to 20)")

MARKER_DICT = {"TUMOR_MARKER": TUMOR_MARKER, "NORMAL_MARKER": NORMAL_MARKER}
max_len = max(len(MARKER_DICT["TUMOR_MARKER"]), len(MARKER_DICT["NORMAL_MARKER"]))
MARKER_DICT["TUMOR_MARKER"] += [None] * (max_len - len(MARKER_DICT["TUMOR_MARKER"]))
MARKER_DICT["NORMAL_MARKER"] += [None] * (max_len - len(MARKER_DICT["NORMAL_MARKER"]))

MARKER_DF = pd.DataFrame(MARKER_DICT)

with pd.ExcelWriter(os.path.join(RESULTS_DIR, "ML_BIOMARKERS_RES.xlsx"), engine="openpyxl") as writer:
    importance_df.to_excel(writer, sheet_name="RANDOM_FOREST_BIOMARKERS", index=False)
    MARKER_DF.to_excel(writer, sheet_name="SHAP_BIOMARKERS", index=False)

tumor_markers_clean = [f for f in MARKER_DICT['TUMOR_MARKER'] if f is not None]
normal_markers_clean = [f for f in MARKER_DICT['NORMAL_MARKER'] if f is not None]

# ==========================================
# 4. SHAP SUMMARY PLOTS SAVING
# ==========================================
selected_feature_indices_1 = [X_test_aligned.columns.get_loc(f) for f in tumor_markers_clean]
X_test_selected_1 = X_test_aligned.iloc[:, selected_feature_indices_1]
shap_values_class_selected_1 = shap_values_class[:, selected_feature_indices_1]

plt.figure(figsize=(10, 8))
shap.summary_plot(
    shap_values_class_selected_1[mask_class_1.values],
    X_test_selected_1[mask_class_1.values],
    plot_type="dot",
    show=False
)
plt.title("Feature Importance - Tumor", fontsize=14)
plt.savefig(os.path.join(FIG_DIR, "TUMOR_FEATURE_IMP_NN_SELECTED.png"), bbox_inches="tight")
plt.close()

selected_feature_indices_0 = [X_test_aligned.columns.get_loc(f) for f in normal_markers_clean]
X_test_selected_0 = X_test_aligned.iloc[:, selected_feature_indices_0]
shap_values_class_selected_0 = shap_values_class[:, selected_feature_indices_0]

plt.figure(figsize=(10, 8))
shap.summary_plot(
    shap_values_class_selected_0[mask_class_0.values],
    X_test_selected_0[mask_class_0.values],
    plot_type="dot",
    show=False
)
plt.title("Feature Importance - Normal", fontsize=14)
plt.savefig(os.path.join(FIG_DIR, "NORMAL_FEATURE_IMP_NN_SELECTED.png"), bbox_inches="tight")
plt.close()

print("Stage 3 complete.")
