# Molecular Property Prediction for Very Large Databases with Natural Language Processing: A Case Study in Ionic Liquid Design

### For calculating the Ionic Liquid properties
There are seceral pretrained ML models for IL properties in the pretrained_models folder and to calculate the properties of ILs, excute the following command
python Mol2Vec_All-IL-Properties-NLP.py

The sample input file is provided in dataset_test folder "filtered_data_IL-LLM_MolecularFormula_Atom-Count.csv"

### For generating novel cations and anions for ILs
Run the python codes Cation-SMILES-Generation_Morpholinium.py.py for generating cation SMILES (example: morpholinium cations) and Anion-SMILES-Generation.py for generating anions SMILES (example: Tf2N anions)

The generated ILs (10.6 millions) and their properties can be found at https://zenodo.org/records/15042249
