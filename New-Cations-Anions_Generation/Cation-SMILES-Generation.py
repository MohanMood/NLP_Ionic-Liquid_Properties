# TF2N-based Anion SMILES Generations

from rdkit import Chem
import csv

def generate_variants(base_smiles, alkyl_chains, functionalities):
    new_molecules = []

    def is_valid_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None

    # Iterate over all combinations of alkyl chains and functionalities
    for chain1 in alkyl_chains:
        try:
            new_smiles_1 = f"[N-](S(=O)(=O){chain1})S(=O)(=O)F"
            if is_valid_smiles(new_smiles_1):
                new_mol_1 = Chem.MolFromSmiles(new_smiles_1)
                new_molecules.append(Chem.MolToSmiles(new_mol_1))
        except Exception as e:
            print(f"Failed to generate molecule with chain {chain1}: {e}")

        for func1 in functionalities:
            try:
                new_smiles_2 = f"[N-](S(=O)(=O){chain1}{func1})S(=O)(=O)F"
                if is_valid_smiles(new_smiles_2):
                    new_mol_2 = Chem.MolFromSmiles(new_smiles_2)
                    new_molecules.append(Chem.MolToSmiles(new_mol_2))
            except Exception as e:
                print(f"Failed to generate molecule with chain {chain1} and functionality {func1}: {e}")

            for chain2 in alkyl_chains:
                try:
                    new_smiles_3 = f"[N-](S(=O)(=O){chain1})S(=O)(=O){chain2}"
                    if is_valid_smiles(new_smiles_3):
                        new_mol_3 = Chem.MolFromSmiles(new_smiles_3)
                        new_molecules.append(Chem.MolToSmiles(new_mol_3))
                except Exception as e:
                    print(f"Failed to generate molecule with chains {chain1}, {chain2}: {e}")

                for func2 in functionalities:
                    try:
                        new_smiles_4 = f"[N-](S(=O)(=O){chain1}{func1})S(=O)(=O){func2}{chain2}"
                        if is_valid_smiles(new_smiles_4):
                            new_mol_4 = Chem.MolFromSmiles(new_smiles_4)
                            new_molecules.append(Chem.MolToSmiles(new_mol_4))

                        new_smiles_5 = f"[N-](S(=O)(=O){func1}{chain1})S(=O)(=O){chain2}{func2}"
                        if is_valid_smiles(new_smiles_5):
                            new_mol_5 = Chem.MolFromSmiles(new_smiles_5)
                            new_molecules.append(Chem.MolToSmiles(new_mol_5))
                    except Exception as e:
                        print(f"Failed to generate molecule with chains {chain1}, {func1}, {chain2}, {func2}: {e}")

        for func1 in functionalities:
            try:
                new_smiles_6 = f"[N-](S(=O)(=O){func1})S(=O)(=O){func2}"
                if is_valid_smiles(new_smiles_6):
                    new_mol_6 = Chem.MolFromSmiles(new_smiles_6)
                    new_molecules.append(Chem.MolToSmiles(new_mol_6))
            except Exception as e:
                print(f"Failed to generate molecule with functionality {func1}: {e}")
    
    return new_molecules

# Base SMILES of Anion (example: NF2)
base_smiles = "[N-](S(=O)(=O)F)S(=O)(=O)F"

# Alkyl chains
alkyl_chains = ["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC", 
    "CCCCCCCCCCC", "CCCCCCCCCCCC", "CCCCCCCCCCCCC", "CCCCCCCCCCCCCC", "CCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCC", 
    "CCCCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCCCC", "C(C)(C)C", "C(C)CC", "C(CC)C", "C(CC)C(C)", "C(C)(CC)C(CC)", "COC", "CCOC", "CCOCC", "CCCOC", "C=C", 
    "C(F)(F)F", "CC(NC(C=C)=O)C", "CSC", "CCSC", "CCCSCC", "CNC", "CCNC", "CCCNC", "CCCNCCC", "CCPC", "CCNCCOC",
    "CCNCCCOCC", "C(CCC)(CC)", "C(CCC)(CCC)CC", "C(C(C)C)(CC)C", "C(C(CC)C)(C(C)C)C", "C(C(CCC)C)(CC(CC)CC)C",
    "C1=CC=CC=C1", "C1=CC=C(N)C=C1", "C1=CC=C(O)C=C1", "C1=CC=C(COC)C=C1", "C1=CC=C(CS(=O)(O)=O)C=C1", "C(F)(F)(F)C(F)(F)C(F)(F)",
    "C1=CC=C(CC(=O)O)C=C1", "C1=CC=C(CC(=O)OC)C=C1", "C(=O)OC", "C1=CC=C(C=C)C=C1", "C1=C(C=C)C=CC=C1", "C(F)(F)(F)C(F)(F)",
    "C1=C(OC)C=CC=C1", "C1=C(NC)C=CC=C1", "CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)", "C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)", "CP(=O)(C)C",
    "CP(=O)(C)OC", "CP(=O)(OC)OC", "CC[Si](OC)(OC)OC", "CC[Si](OC)(OC)C", "CC[Si](OC)(C)C", "COCCOC", "C(F)(F)(F)", 
   "O=C(C)OC", "OCCN(CO)CC", "N#C", "CC(NC(C=C)=O)C", "C(N(CCO)CCO)C", "N#CC", "O[C@H](CN1CCOCC1)C", "C(N(CCO)CCO)C",
   "C[C@@H](CC1=CC=CC=C1)N", "C[C@@H](CC1=CC=CC=C1)NC", "CC1=CC=C(C#N)C=C1", "CC1=CC(C#N)=CC=C1", "CC1=CC=C(F)C=C1",
   "CC1=CC(F)=CC=C1", "CC1=C(F)C=CC=C1", "CC1=C(C(F)(F)F)C=CC=C1", "CC1=CC=C(C(F)(F)F)C=C1", "CC1=CC(OC)=CC=C1", "Cc1ccc(F)c1"
   "CC1=C(OC)C=CC=C1", "CC1=CC=C(OC)C=C1", "CC1=CC=CC=C1", "C1=CC2=C(C=CC=C2)C=C1", "C[C@@H](CC1=CC=C(N)C=C1)NC", "CCCCCCCC/C=C\CCCCCCCC",
   "CN1C=CN=C1", "C1=NC=CN1C", "C1N=CN(C)C=1", "C1N(C=NC=1)C", "N1=NC=CN1", "C1NN=NC=1", "N1N=NN=C1", "N1(CC(C)C)N=NN=C1", "N1=NN(C2=CC=CC=C2)N=C1", 
   "N1=NN(C2=CC(OC)=C(N)C(OC)=C2)N=C1", "N1=NN(C2=CC(OC)=C(OC)C(OC)=C2)N=C1", "N1=NN(C2=CC(OC)=C(OC)C(OC)=C2)N=C1", "C1=CC(OC)=C(OC)C(OC)=C1", 
    "C1=CC(OC)=C(O)C(OC)=C1", "CN1C(C2=CC=CC=C2)=NC=C1", "CN1C(C2=CC(OC)=CC(OC)=C2)=NC=C1", "CN1C(C2=CC(O)=CC(O)=C2)=NC=C1", "C1=CC=C(O)C(OC)=C1", "C1=CC=C(O)C=C1",
    "CC(=O)[NH]C1C=COC=1", "O1C=C([NH]C(C)=O)C=C1", "O1C=CC=C1", "C1=COC=C1", "COC(C1=CC=CC=C1)", "COC(C1=CC(OC)=CC=C1)", "COC(C1=CC(OC)=CC(OC)=C1)", "COC(C1=CC(OC)=C(O)C(OC)=C1)",
    "COC(C1=CC(OC)=C(N)C(OC)=C1)", "CNC(C1=CC(OC)=C(N)C(OC)=C1)", "CNC(C1=CC(OC)=C(O)C(OC)=C1)", "C[N+]([O-])=O", "C[N+](C)(O)C", "S(=O)(=O)O", "S(=O)(=O)(O)C"]

# Functionalities
functionalities = ["O", "N", "Cl", "Br", "S"]

# Generate new molecules
new_molecules = generate_variants(base_smiles, alkyl_chains, functionalities)

# Print new SMILES strings
for smiles in new_molecules:
    print(smiles)

# Save new SMILES strings to a CSV file
with open('TF2N-based_Anions-SMILES_2.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["SMILES"])
    for smiles in new_molecules:
        writer.writerow([smiles])
