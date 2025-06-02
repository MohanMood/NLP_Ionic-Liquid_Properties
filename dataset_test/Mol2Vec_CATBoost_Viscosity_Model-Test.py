import sys
import os
sys.path.append('path_of_iSIM/iSIM_Clustering')

import math
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import seaborn as sns
from mol2vec.features import mol2alt_sentence, MolSentence, DfVec, sentences2vec
from mol2vec.helpers import depict_identifier, plot_2D_vectors, IdentifierTable, mol_to_svg
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.model_selection import KFold
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from xgboost import XGBRegressor
from gensim.models import word2vec
from sklearn import preprocessing
from catboost import CatBoostRegressor
import optuna
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
from div import *
from sampling import *
from comp import *
from real import calculate_isim_real, calculate_medoid_real, calculate_outlier_real, pairwise_average_real
from utils import binary_fps, minmax_norm, real_fps, pairwise_average, rdkit_pairwise_sim
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from IPython.display import display, HTML

# Reading Ionic Liquids Viscosity Data
dataset_csv = 'IL-Viscosity'
input_dir = 'path_of_dataset/Dataset/'
dataset = pd.read_csv(input_dir + f'{dataset_csv}.csv', encoding= 'unicode_escape')
print("Viscosity Dataset:", dataset.shape)

smiles = dataset['IL_SMILES'].tolist()

# Generation of fingerprints
fps = binary_fps(smiles, fp_type='ECFP4', n_bits=256)

# Optional: save the fingerprints in a npy file
fp_output_path = os.path.join(input_dir, 'ILs_Viscosity.npy')
np.save(fp_output_path, fps)

# === Load fingerprints and define output for t-SNE ===
fingerprints = np.load(fp_output_path)
tSNE_output = os.path.join('tSNE_samp.csv')

# Define the similarity index for sampling
n_ary = 'JT'

# Define the percentage for initial data split (80% training and 20% testing)
percentage = 80

# Get the sampled indexes
strat_samp = stratified_sampling(fingerprints, n_ary, percentage)

# Do tSNE to visualize the sampled indexes
tsne = TSNE(n_components=2)
X =tsne.fit_transform(fingerprints)

# Create a DataFrame with tSNE coordinates
tsne_df = pd.DataFrame(X, columns=['tSNE1', 'tSNE2'])

label_1 = 'Testing'
label_2 = 'Training (Sampling)'
label_3 = 'Validation'

# Initialize the label column as 'testing'
tsne_df['Label'] = label_1

# Label Sampling samples (80% of data)
tsne_df.loc[strat_samp, 'Label'] = label_2

# Extract Sampling subset (80% of data)
dataset = pd.read_csv(input_dir + f'{dataset_csv}.csv', encoding='unicode_escape')
dataset = pd.concat([dataset, tsne_df], axis=1)
dataset_sampling = dataset[dataset['Label'] == label_2].copy()

# Generate fingerprints for Sampling subset
smiles_sampling = dataset_sampling['IL_SMILES'].tolist()
fps_sampling = binary_fps(smiles_sampling, fp_type='ECFP4', n_bits=256)
np.save(input_dir + 'IL_Viscosity_Sampling.npy', fps_sampling)
fingerprints_sampling = np.load(input_dir + 'IL_Viscosity_Sampling.npy')

# Further split data into training (70%) and validation (10%)
percentage_new = 87.5  # 70/80 = 0.875
sampling_new = stratified_sampling(fingerprints_sampling, n_ary, percentage_new)

# Create training and validation indices
train_indices = dataset_sampling.index[sampling_new].tolist()
val_indices = dataset_sampling.index.difference(train_indices).tolist()

# Update tsne_df['Label'] to reflect training and validation splits
tsne_df.loc[train_indices, 'Label'] = label_2
tsne_df.loc[val_indices, 'Label'] = label_3
tsne_df.loc[dataset[dataset['Label'] == label_1].index, 'Label'] = label_1

# Save the updated tSNE DataFrame
tsne_df.to_csv(input_dir + tSNE_output, index=False)
tsne_df = pd.read_csv(input_dir + tSNE_output, encoding='unicode_escape')

# Combine dataset with updated tsne_df
dataset = pd.concat([dataset.drop(['tSNE1', 'tSNE2', 'Label'], axis=1, errors='ignore'), tsne_df], axis=1)


# Extracting necessary columns
Smiles_IL = dataset['IL_SMILES']
Smiles_IL_Unique = dataset['IL_SMILES'].unique()
print("Unique IL SMILES:", len(Smiles_IL_Unique))

target_prop = dataset['ln_visc']
temperature_prop = dataset['Temp (K)']
pressure_prop = dataset['Pressure (kPa)']


# Load the pre-trained Word2Vec model
model_path = 'path_of_mol2vec/model_300dim.pkl'             # Install Mol2vec (https://github.com/samoturk/mol2vec)
model = word2vec.Word2Vec.load(model_path)

# Extraction of identifiers from molecules
IL_smiles_prop = [Chem.MolFromSmiles(x) for x in Smiles_IL]
IL_sentences_prop = [mol2alt_sentence(x, 1) for x in IL_smiles_prop]

# Define the DfVec class if not already defined
class DfVec:
    def __init__(self, vector):
        self.vector = vector

# Convert sentences to vectors
def sentences2vec(sentences, model, unseen='UNK'):
    vectors = []
    for sentence in sentences:
        vec = []
        for word in sentence:
            if word in model.wv.key_to_index:
                vec.append(model.wv[word])
            else:
                vec.append(model.wv[unseen] if unseen in model.wv.key_to_index else np.zeros(model.vector_size))
        vectors.append(np.sum(vec, axis=0))
    return vectors

# Convert sentences to vectors
IL_mol2vec_prop = [DfVec(x) for x in sentences2vec(IL_sentences_prop, model, unseen='UNK')]
np_IL_mol2vec_prop = np.array([x.vector for x in IL_mol2vec_prop])

# Standardize the mol2vec features
scaler = preprocessing.StandardScaler()
np_IL_mol2vec_prop_scaler = scaler.fit_transform(np_IL_mol2vec_prop)

# Concatenate atom counts to the mol2vec features
np_IL_features_prop = np.concatenate((np_IL_mol2vec_prop, temperature_prop.values.reshape(-1, 1), 
                                      pressure_prop.values.reshape(-1, 1)), axis=1)

print("features shape:", np_IL_features_prop.shape)


# Assign training, validation, and test sets
X_train = np_IL_features_prop[train_indices]
y_train = target_prop[train_indices]
smiles_train = Smiles_IL[train_indices]
temp_train = temperature_prop[train_indices]
press_train = pressure_prop[train_indices]

X_val = np_IL_features_prop[val_indices]
y_val = target_prop[val_indices]
smiles_val = Smiles_IL[val_indices]
temp_val = temperature_prop[val_indices]
press_val = pressure_prop[val_indices]

X_test = np_IL_features_prop[dataset[dataset['Label'] == label_1].index]
y_test = target_prop[dataset[dataset['Label'] == label_1].index]
smiles_test = Smiles_IL[dataset[dataset['Label'] == label_1].index]
temp_test = temperature_prop[dataset[dataset['Label'] == label_1].index]
press_test = pressure_prop[dataset[dataset['Label'] == label_1].index]



# Train CatBoost with optimized hyperparameters
cat_model = CatBoostRegressor(verbose=100, n_estimators = 2798, l2_leaf_reg = 6, learning_rate = 0.17741819704159587, 
                              depth = 5, random_strength = 1, rsm = 0.7875510288973779, boosting_type='Plain', 
                              early_stopping_rounds=30)

cat_model.fit(X_train, y_train, eval_set=(X_val, y_val))

# Predict on train
y_pred_train = cat_model.predict(X_train)

# Print metrics
print("\nTrain MAE:", mean_absolute_error(y_train, y_pred_train))
print("Train RMSE:", math.sqrt(mean_squared_error(y_train, y_pred_train)))
print("Train R2:", r2_score(y_train, y_pred_train))

# Predict on Validation
y_pred_val = cat_model.predict(X_val)

# Print metrics
print("\nValidation MAE:", mean_absolute_error(y_val, y_pred_val))
print("Validation RMSE:", math.sqrt(mean_squared_error(y_val, y_pred_val)))
print("Validation R2:", r2_score(y_val, y_pred_val))

# Predict on test
y_pred_test = cat_model.predict(X_test)

# Print metrics
print("\nTest MAE:", mean_absolute_error(y_test, y_pred_test))
print("Test RMSE:", math.sqrt(mean_squared_error(y_test, y_pred_test)))
print("Test R2:", r2_score(y_test, y_pred_test))

# Create separate DataFrames for each split
train_df = pd.DataFrame({
    'SMILES': smiles_train,
    'Temp (K)': temp_train,
    'viscosity_exp': y_train,
    'viscosity_pred': y_pred_train,
    'Set': 'Training'
})

val_df = pd.DataFrame({
    'SMILES': smiles_val,
    'Temp (K)': temp_val,
    'viscosity_exp': y_val,
    'viscosity_pred': y_pred_val,
    'Set': 'Validation'
})

test_df = pd.DataFrame({
    'SMILES': smiles_test,
    'Temp (K)': temp_test,
    'viscosity_exp': y_test,
    'viscosity_pred': y_pred_test,
    'Set': 'Testing'
})

# Combine them into one
combined_df = pd.concat([train_df, val_df, test_df], ignore_index=True)

# Save to CSV
combined_df.to_csv('CATBoost_Mol2vec-Viscosity_pred.csv', index=False)

# Save the trained model to a file
cat_model.save_model('CATBoost_Mol2vec_Viscosity_Model_Stratified.cbm')
print("Trained Viscosity Model Saved Succesfully")
