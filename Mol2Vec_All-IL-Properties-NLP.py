import math
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
from mol2vec.features import mol2alt_sentence
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.model_selection import KFold
from sklearn import preprocessing
from catboost import CatBoostRegressor, CatBoostClassifier
from gensim.models import word2vec

# Pretrained Models Path
pretrained_models = 'Pretrained_model_path/'

# Load the saved models
model_viscosity = CatBoostRegressor()
model_viscosity.load_model(pretrained_models + 'CATBoost_Mol2vec_Viscosity_Model_Stratified.cbm')

model_conductivity = CatBoostRegressor()
model_conductivity.load_model(pretrained_models + 'CATBoost_Mol2vec_Ionic_Conductivity_model_Stratified.cbm')

model_density = CatBoostRegressor()
model_density.load_model(pretrained_models + 'CATBoost_Mol2vec-Density_model_Stratified.cbm')

model_surface = CatBoostRegressor()
model_surface.load_model(pretrained_models + 'CATBoost_Mol2vec-ST-Stratified.cbm')

model_toxicity = CatBoostRegressor()
model_toxicity.load_model(pretrained_models + 'CATBoost_Mol2vec-toxicity-Stratified.cbm')

model_melting = CatBoostClassifier()
model_melting.load_model(pretrained_models + 'CATBoost_Mol2vec_Melting_Stratified_Classification.cbm')

model_activity_class = CatBoostClassifier()
model_activity_class.load_model(pretrained_models + 'CATBoost_Mol2vec_Water_Activity_Stratified_Classification.cbm')

model_activity_reg = CatBoostRegressor()
model_activity_reg.load_model(pretrained_models + 'CATBoost_Mol2vec-Water_Activity_model_Stratified.cbm')


# Predictions on the Novel IL Properties without further Training a Models
# New dataset for predictions
input_dir_new = 'input_directory_path/'
dataset_new = pd.read_csv(input_dir_new + 'filtered_data_IL-LLM_MolecularFormula_Atom-Count.csv', encoding='unicode_escape')

# Extracting necessary columns for new dataset
IL_smiles = dataset_new['SMILES']
temperature_new = dataset_new['Temperature, K']
pressure_new = dataset_new['Pressure, kPa']
water_molfrac_activity_new = dataset_new['Mole_frac_water']
atom_counts_new = dataset_new[['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Si']].values

# Load the pre-trained Word2Vec model (renamed variable to avoid conflict with CatBoostRegressor)
word2vec_model_path = 'pre-trained_mol2vec_path/model_300dim.pkl'
model = word2vec.Word2Vec.load(word2vec_model_path)

# Extraction of identifiers from molecules
IL_smiles_new = [Chem.MolFromSmiles(x) for x in IL_smiles]
IL_sentences_new = [mol2alt_sentence(x, 1) for x in IL_smiles_new]

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
IL_mol2vec_new = [DfVec(x) for x in sentences2vec(IL_sentences_new, model, unseen='UNK')]
np_IL_mol2vec_new = np.array([x.vector for x in IL_mol2vec_new])

# Standardize the mol2vec features
scaler = preprocessing.StandardScaler()
np_IL_mol2vec_new_scaler = scaler.fit_transform(np_IL_mol2vec_new)

# Concatenate atom counts to the new mol2vec features
input_features_new = np.concatenate((np_IL_mol2vec_new, temperature_new.values.reshape(-1, 1), pressure_new.values.reshape(-1, 1)), axis=1)
input_features_surface_new = np.concatenate((np_IL_mol2vec_new, temperature_new.values.reshape(-1, 1)), axis=1)
input_features_new_activity = np.concatenate((np_IL_mol2vec_new, 
                                     water_molfrac_activity_new.values.reshape(-1, 1),
                                     temperature_new.values.reshape(-1, 1), 
                                     pressure_new.values.reshape(-1, 1)), axis=1)

# Predict viscosity for the new dataset using the already trained CatBoost model
target_pred_visc_new = model_viscosity.predict(input_features_new)
target_pred_cond_new = model_conductivity.predict(input_features_new)
target_pred_density_new = model_density.predict(input_features_new)
target_pred_surface_new = model_surface.predict(input_features_surface_new)
target_pred_toxicity_new = model_toxicity.predict(np_IL_mol2vec_new)
target_pred_melting_new = model_melting.predict(np_IL_mol2vec_new)
target_pred_activity_class = model_activity_class.predict(input_features_new_activity)

# Flatten the atom counts for each molecule and save the testing results
atom_test_flat = [','.join(map(str, at)) for at in atom_counts_new]
results_df = pd.DataFrame({
    'Canonical_SMILES': IL_smiles_new,
    'Temp (K)': temperature_new,
    'Pressure (kPa)': pressure_new,
    'mole_fraction': water_molfrac_activity_new,
    'ln(viscosity), mPa.s': target_pred_visc_new,
    'Density, kg/m3': target_pred_density_new,
    'Ionic conductivity, S/m': target_pred_cond_new,
    'Surface tension, mN/m': target_pred_surface_new,
    'Toxicity (logEC50)': target_pred_toxicity_new,
    'Melting Temp. (K)': target_pred_melting_new,
    'Water activity in ILs': target_pred_activity_class
})

Output_dir = 'output_directory_path/'
results_df.to_csv(Output_dir + f'Predictions_filtered_data.csv', index=False)

print(" ")
print("Predictions Done and the Data Was Saved Succesfully")
print(" ")
