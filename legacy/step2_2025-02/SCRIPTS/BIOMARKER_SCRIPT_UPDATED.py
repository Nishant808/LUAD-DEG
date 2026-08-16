# Import necessary libraries
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
from sklearn.preprocessing import StandardScaler
import os
import shap
import numpy as np
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
import matplotlib
matplotlib.use('Agg')  # Use 'Agg' backend for non-interactive plotting
import matplotlib.pyplot as plt

# set Directories

MAIN_DIR = "Cancer_Biomarkers_LUAD/ML_Implemenation_Step-2"

try:
    os.chdir(MAIN_DIR)
except:
    print("SKIP") 

plot_dir = os.path.join(MAIN_DIR,"BIOMARKERS_RES")

try:
    os.mkdir(plot_dir)
except:
    print("EXIST")


    
# Load the data (assuming your data is in a CSV file)
data = pd.read_csv('DATA/LAUD_FEATURE_DATA_MATRIX.csv', index_col=0)

# Split data into features and target variable
X = data.drop('Label', axis=1)
y = data['Label']

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Convert the scaled features back to a DataFrame
X_scaled_df = pd.DataFrame(X_scaled, columns=X.columns)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled_df, y, test_size=0.2, random_state=42)


# Define the model
model = Sequential([
    Dense(128, input_dim=X_train.shape[1], activation='relu'),
    Dropout(0.6),  # Dropout for regularization
    Dense(64, activation='relu'),  # Hidden layer
    Dropout(0.6),  # Dropout for regularization
    Dense(32, activation='relu'),  # Hidden layer
    Dropout(0.6),  # Dropout for regularization
    Dense(16, activation='relu'),  # Hidden layer
    Dropout(0.6),  # Dropout for regularization
    Dense(1, activation='sigmoid')  # Output layer for binary classification
])
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.fit(X_train, y_train, epochs=100, batch_size=32, validation_split=0.2)

X_train_np = X_train.values if isinstance(X_train, pd.DataFrame) else np.array(X_train)
X_test_np = X_test.values if isinstance(X_test, pd.DataFrame) else np.array(X_test)

explainer = shap.DeepExplainer(model, X_train_np)

# Calculate SHAP values for the test set
shap_values = explainer.shap_values(X_test_np)

# For binary classification, shap_values[0] is used
shap_values_class = shap_values[0]

# Ensure indices of y_test and X_test are aligned
y_test_aligned = y_test.reset_index(drop=True)
X_test_aligned = X_test.reset_index(drop=True)

mask_class_0 = y_test_aligned == 0
mask_class_1 = y_test_aligned == 1

shap_values_class_1 = shap_values_class[mask_class_1.values]
feature_importances_class_1 = np.abs(shap_values_class_1).mean(0)
feature_names_class_1 = X_test_aligned.columns.tolist()

SHAP_CLASS_1 = pd.DataFrame({
    'Features': feature_names_class_1,
    'SHAP Scores': feature_importances_class_1
})



shap_values_class_0 = shap_values_class[mask_class_0.values]
feature_importances_class_0 = np.abs(shap_values_class_0).mean(0)
feature_names_class_0 = X_test_aligned.columns.tolist()
SHAP_CLASS_0 = pd.DataFrame({
    'Features': feature_names_class_0,
    'SHAP Scores': feature_importances_class_0
})

# Sort the DataFrame by SHAP scores in descending order
SHAP_CLASS_0 = SHAP_CLASS_0.sort_values(by='SHAP Scores', ascending=False)
SHAP_CLASS_1 = SHAP_CLASS_1.sort_values(by='SHAP Scores', ascending=False)


MARKER_DICT = {
    "TUMOR_MARKER" : [],
    "NORMAL_MARKER" : []
}

for G_1, G_0 in zip(SHAP_CLASS_1['Features'], SHAP_CLASS_0['Features']):
    if (len(MARKER_DICT["TUMOR_MARKER"]) >= 20) and (len(MARKER_DICT["NORMAL_MARKER"]) >= 20):
        break
    if G_1 not in MARKER_DICT["NORMAL_MARKER"]:
        MARKER_DICT["TUMOR_MARKER"].append(G_1)
    else:
        MARKER_DICT["NORMAL_MARKER"].remove(G_1)
        
    if G_0 not in MARKER_DICT["TUMOR_MARKER"]:
        MARKER_DICT["NORMAL_MARKER"].append(G_0)
    else:
        MARKER_DICT["TUMOR_MARKER"].remove(G_0)
        

MARKER_DF = pd.DataFrame(MARKER_DICT)
with pd.ExcelWriter("BIOMARKERS_RES/ML_BIOMARKERS_RES.xlsx","openpyxl") as writter:
    MARKER_DF.to_excel(writter,sheet_name="SHAP_BIOMARKERS", index=False)


# Filter the SHAP values and X_test for the selected features
selected_feature_indices_1 = [X_test_aligned.columns.get_loc(f) for f in MARKER_DICT['TUMOR_MARKER']]
X_test_selected_1 = X_test_aligned.iloc[:, selected_feature_indices_1]
shap_values_class_selected_1 = shap_values_class[:, selected_feature_indices_1]

# Plot SHAP summary plot for Class 1 (assuming Tumor class)
plt.figure(figsize=(10, 8))
shap.summary_plot(shap_values_class_selected_1[mask_class_1.values], X_test_selected_1[mask_class_1.values], 
                  plot_type="dot", title="Feature Importance - Tumor")
plt.savefig("BIOMARKERS_RES/TUMOR_FEATURE_IMP_NN_SELECTED.png", bbox_inches="tight")
plt.close()  # Close the figure to ensure it doesn't interfere with the next plot


# Filter the SHAP values and X_test for the selected features
selected_feature_indices_0 = [X_test_aligned.columns.get_loc(f) for f in MARKER_DICT["NORMAL_MARKER"]]
X_test_selected_0 = X_test_aligned.iloc[:, selected_feature_indices_0]
shap_values_class_selected_0 = shap_values_class[:, selected_feature_indices_0]


# Plot SHAP summary plot for Class 0 (assuming Normal class)
plt.figure(figsize=(10, 8))
shap.summary_plot(shap_values_class_selected_0[mask_class_0.values], X_test_selected_0[mask_class_0.values], 
                  plot_type="dot", title="Feature Importance - Normal")
plt.savefig("BIOMARKERS_RES/NORMAL_FEATURE_IMP_NN_SELECTED.png", bbox_inches="tight")
plt.close()  # Close the figure to ensure it doesn't interfere with the next plot

