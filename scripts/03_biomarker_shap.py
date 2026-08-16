# Import necessary libraries
import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')  # Use 'Agg' backend for non-interactive plotting
import matplotlib.pyplot as plt
import shap

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
from sklearn.preprocessing import StandardScaler

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout



# Set Directories
MAIN_DIR = "/Users/nishantthalwal/Documents/Data ALL/TCGA_LUAD"

try:
    os.chdir(MAIN_DIR)
except Exception as e:
    print(f"SKIP directory change: {e}")

plot_dir = os.path.join(MAIN_DIR, "results", "biomarkers")

try:
    os.makedirs(plot_dir, exist_ok=True)
except Exception as e:
    print(f"Directory warning: {e}")

# Load the data
data = pd.read_csv(os.path.join(MAIN_DIR, 'data', 'FEATURE_DATA_MATRIX.csv'), index_col=0)
print("Data shape:", data.shape)

# Split data into features and target variable
X = data.drop('Label', axis=1)
y = data['Label']

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Convert the scaled features back to a DataFrame
X_scaled_df = pd.DataFrame(X_scaled, columns=X.columns)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled_df, y, test_size=0.2, random_state=42)

# ==========================================
# 1. RANDOM FOREST MODEL & FEATURE IMPORTANCE
# ==========================================
rf_model = RandomForestClassifier(n_estimators=100, random_state=52)
rf_model.fit(X_train, y_train)

# Make predictions
y_pred = rf_model.predict(X_test)

# Extract feature importances
feature_importances = rf_model.feature_importances_
feature_names = X.columns

# Create a DataFrame for visualization
importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': feature_importances})
importance_df = importance_df.sort_values(by='Importance', ascending=False)

# Plot feature importances
plt.figure(figsize=(12, 8))
plt.barh(importance_df['Feature'][:20], importance_df['Importance'][:20])
plt.xlabel('Importance')
plt.ylabel('Feature')
plt.title('Top 20 Genes : RandomForest')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.gca().invert_yaxis()  # Invert y-axis to show top features on top
plt.savefig("results/biomarkers/Top_MARKERS_RandomForest.png", bbox_inches="tight")
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
model.fit(X_train, y_train, epochs=50, batch_size=32, validation_split=0.2, verbose=1)

X_train_np = X_train.values if isinstance(X_train, pd.DataFrame) else np.array(X_train)
X_test_np = X_test.values if isinstance(X_test, pd.DataFrame) else np.array(X_test)

# Use generic shap.Explainer to handle modern Keras/TF execution styles natively
explainer = shap.Explainer(model, X_train_np)

# FIX: Added max_evals=1000 to allow the Permutation Explainer enough budget for 309 features
shap_values = explainer(X_test_np, max_evals=1000)

# Extract the raw output values array
if hasattr(shap_values, "values"):
    shap_values_matrix = shap_values.values
else:
    shap_values_matrix = np.array(shap_values)

# Fix dimensional shifts (ensures shape matches (samples, features))
if len(shap_values_matrix.shape) == 3:
    shap_values_matrix = np.squeeze(shap_values_matrix, axis=-1)
elif shap_values_matrix.shape[0] != X_test_np.shape[0] and shap_values_matrix.shape[1] == X_test_np.shape[0]:
    shap_values_matrix = shap_values_matrix.T

shap_values_class = shap_values_matrix

# Ensure indices of y_test and X_test are aligned
y_test_aligned = y_test.reset_index(drop=True)
X_test_aligned = X_test.reset_index(drop=True)

mask_class_0 = (y_test_aligned == 0)
mask_class_1 = (y_test_aligned == 1)

# Process Class 1 (Tumor)
shap_values_class_1 = shap_values_class[mask_class_1.values]
feature_importances_class_1 = np.abs(shap_values_class_1).mean(0)
SHAP_CLASS_1 = pd.DataFrame({
    'Features': X_test_aligned.columns.tolist(),
    'SHAP Scores': feature_importances_class_1
}).sort_values(by='SHAP Scores', ascending=False)

# Process Class 0 (Normal)
shap_values_class_0 = shap_values_class[mask_class_0.values]
feature_importances_class_0 = np.abs(shap_values_class_0).mean(0)
SHAP_CLASS_0 = pd.DataFrame({
    'Features': X_test_aligned.columns.tolist(),
    'SHAP Scores': feature_importances_class_0
}).sort_values(by='SHAP Scores', ascending=False)

# ==========================================
# 3. BIOMARKER SELECTION LOGIC
# ==========================================
MARKER_DICT = {
    "TUMOR_MARKER": [],
    "NORMAL_MARKER": []
}

for G_1, G_0 in zip(SHAP_CLASS_1['Features'], SHAP_CLASS_0['Features']):
    if len(MARKER_DICT["TUMOR_MARKER"]) < 20:
        if G_1 not in MARKER_DICT["NORMAL_MARKER"]:
            MARKER_DICT["TUMOR_MARKER"].append(G_1)
        elif G_1 in MARKER_DICT["NORMAL_MARKER"]:
            MARKER_DICT["NORMAL_MARKER"].remove(G_1)

    if len(MARKER_DICT["NORMAL_MARKER"]) < 20:
        if G_0 not in MARKER_DICT["TUMOR_MARKER"]:
            MARKER_DICT["NORMAL_MARKER"].append(G_0)
        elif G_0 in MARKER_DICT["TUMOR_MARKER"]:
            MARKER_DICT["TUMOR_MARKER"].remove(G_0)

    if len(MARKER_DICT["TUMOR_MARKER"]) >= 20 and len(MARKER_DICT["NORMAL_MARKER"]) >= 20:
        break

# Pad uneven lengths to prevent DataFrame crash
max_len = max(len(MARKER_DICT["TUMOR_MARKER"]), len(MARKER_DICT["NORMAL_MARKER"]))
MARKER_DICT["TUMOR_MARKER"] += [None] * (max_len - len(MARKER_DICT["TUMOR_MARKER"]))
MARKER_DICT["NORMAL_MARKER"] += [None] * (max_len - len(MARKER_DICT["NORMAL_MARKER"]))

MARKER_DF = pd.DataFrame(MARKER_DICT)

# Save to Excel
with pd.ExcelWriter("results/biomarkers/ML_BIOMARKERS_RES.xlsx", engine="openpyxl") as writer:
    importance_df.to_excel(writer, sheet_name="RANDOM_FOREST_BIOMARKERS", index=False)
    MARKER_DF.to_excel(writer, sheet_name="SHAP_BIOMARKERS", index=False)

# Filter out padding None values for reliable plotting mapping
tumor_markers_clean = [f for f in MARKER_DICT['TUMOR_MARKER'] if f is not None]
normal_markers_clean = [f for f in MARKER_DICT['NORMAL_MARKER'] if f is not None]

# ==========================================
# 4. SHAP SUMMARY PLOTS SAVING
# ==========================================
# Plot Tumor Class Summary
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
plt.savefig("results/biomarkers/TUMOR_FEATURE_IMP_NN_SELECTED.png", bbox_inches="tight")
plt.close()

# Plot Normal Class Summary
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
plt.savefig("results/biomarkers/NORMAL_FEATURE_IMP_NN_SELECTED.png", bbox_inches="tight")
plt.close()