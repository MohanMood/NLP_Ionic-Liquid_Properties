from rdkit import Chem
import csv

def generate_variants(smiles, alkyl_chains, functionalities):
    new_molecules = []

    def is_valid_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None

    # Iterate over all combinations of alkyl chains at specified positions
    for chain1 in alkyl_chains:
        for chain2 in alkyl_chains:
            for chain3 in alkyl_chains:
                for chain4 in alkyl_chains:
                    try:
                        # Modify SMILES to include the alkyl chains at specified positions
                        new_smiles = f"{chain1}C[N+]1(C{chain2})CC(C{chain3})OC(C{chain4})C1"
                        if is_valid_smiles(new_smiles):
                            new_mol = Chem.MolFromSmiles(new_smiles)
                            new_molecules.append(Chem.MolToSmiles(new_mol))

                        # Add alkyl chains and functionalities at the specified positions
                        for func in functionalities:
                            new_smiles_func_1 = f"{chain1}C[N+]1(C)CCOCC1"
                            if is_valid_smiles(new_smiles_func_1):
                                new_mol_func_1 = Chem.MolFromSmiles(new_smiles_func_1)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_1))

                            new_smiles_func_2 = f"{func}{chain1}C[N+]1(C)CCOCC1"
                            if is_valid_smiles(new_smiles_func_2):
                                new_mol_func_2 = Chem.MolFromSmiles(new_smiles_func_2)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_2))

                            new_smiles_func_3 = f"{func}{chain1}C[N+]1(C{chain2})CCOCC1"
                            if is_valid_smiles(new_smiles_func_3):
                                new_mol_func_3 = Chem.MolFromSmiles(new_smiles_func_3)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_3))                                    

                            new_smiles_func_4 = f"{func}{chain1}C[N+]1(C{func}{chain2})CCOCC1"
                            if is_valid_smiles(new_smiles_func_4):
                                new_mol_func_4 = Chem.MolFromSmiles(new_smiles_func_4)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_4))                                    

                            new_smiles_func_5 = f"{func}{chain1}C[N+]1(C{func}{chain2}{func})CCOCC1"
                            if is_valid_smiles(new_smiles_func_5):
                                new_mol_func_5 = Chem.MolFromSmiles(new_smiles_func_5)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_5))

                            new_smiles_func_6 = f"{chain1}C[N+]1(C{chain2})CCOCC1"
                            if is_valid_smiles(new_smiles_func_6):
                                new_mol_func_6 = Chem.MolFromSmiles(new_smiles_func_6)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_6))

                            new_smiles_func_7 = f"{chain1}C[N+]1(C{chain2})CC(C{chain3})OCC1"
                            if is_valid_smiles(new_smiles_func_7):
                                new_mol_func_7 = Chem.MolFromSmiles(new_smiles_func_7)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_7))

                            new_smiles_func_8 = f"{func}{chain1}C[N+]1(C{chain2})CC(C{chain3})OCC1"
                            if is_valid_smiles(new_smiles_func_8):
                                new_mol_func_8 = Chem.MolFromSmiles(new_smiles_func_8)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_8))

                            new_smiles_func_9 = f"{func}{chain1}C[N+]1(C{func}{chain2})CC(C{chain3})OCC1"
                            if is_valid_smiles(new_smiles_func_9):
                                new_mol_func_9 = Chem.MolFromSmiles(new_smiles_func_9)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_9))

                            new_smiles_func_10 = f"{func}{chain1}C[N+]1(C{func}{chain2})CC(C{func}{chain3})OCC1"
                            if is_valid_smiles(new_smiles_func_10):
                                new_mol_func_10 = Chem.MolFromSmiles(new_smiles_func_10)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_10))

                            new_smiles_func_11 = f"{func}{chain1}{func}C[N+]1(C{chain2}{func})CC(C{func}{chain3})OCC1"
                            if is_valid_smiles(new_smiles_func_11):
                                new_mol_func_11 = Chem.MolFromSmiles(new_smiles_func_11)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_11))

                            new_smiles_func_12 = f"{func}{chain1}{func}C[N+]1(C{chain2}{func})CC(C{chain3})OCC1"
                            if is_valid_smiles(new_smiles_func_12):
                                new_mol_func_12 = Chem.MolFromSmiles(new_smiles_func_12)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_12))

                            new_smiles_func_13 = f"{chain1}{func}C[N+]1(C{chain2})CC(C{chain3})OC(C{chain4})C1"
                            if is_valid_smiles(new_smiles_func_13):
                                new_mol_func_13 = Chem.MolFromSmiles(new_smiles_func_13)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_13))

                            new_smiles_func_14 = f"{chain1}{func}C[N+]1(C{func}{chain2})CC(C{func}{chain3})OC(C{chain4})C1"
                            if is_valid_smiles(new_smiles_func_14):
                                new_mol_func_14 = Chem.MolFromSmiles(new_smiles_func_14)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_14))

                            new_smiles_func_15 = f"{chain1}{func}C[N+]1(C{func}{chain2})CC(C{func}{chain3})OC(C{func}{chain4})C1"
                            if is_valid_smiles(new_smiles_func_15):
                                new_mol_func_15 = Chem.MolFromSmiles(new_smiles_func_15)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_15))

                            new_smiles_func_16 = f"{chain1}{func}C[N+]1(C{func}{chain2})CC(C{func}{chain3})OC(C{chain4}{func})C1"
                            if is_valid_smiles(new_smiles_func_16):
                                new_mol_func_16 = Chem.MolFromSmiles(new_smiles_func_16)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_16))

                            new_smiles_func_17 = f"{chain1}{func}C[N+]1(C{chain2})CC(C{chain3})OC(C{chain4}{func})C1"
                            if is_valid_smiles(new_smiles_func_17):
                                new_mol_func_17 = Chem.MolFromSmiles(new_smiles_func_17)
                                new_molecules.append(Chem.MolToSmiles(new_mol_func_17))

                    except Exception as e:
                        print(f"Failed to generate molecule with chains {chain1}, {chain2}, {chain4}, and {chain5}: {e}")

    return new_molecules

# Base SMILES of Cation
base_smiles = "C[N+]1(C)CCOCC1"

# Alkyl chains
alkyl_chains = ["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC", 
    "CCCCCCCCCCC", "CCCCCCCCCCCC", "CCCCCCCCCCCCC", "CCCCCCCCCCCCCC", "CCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCC", 
    "CCCCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCCCC", "C(C)(C)C", "C(C)CC", "C(CC)C", "C(CC)C(C)", "C(C)(CC)C(CC)", "COC", "CCOC", "CCOCC", "CCCOC", "C=C", 
    "C(F)(F)F", "CC(NC(C=C)=O)C", "CSC", "CCSC", "CCCSCC", "CNC", "CCNC", "CCCNC", "CCCNCCC", "CCPC", "CCNCCOC",
    "CCNCCCOCC", "C(CCC)(CC)", "C(CCC)(CCC)CC", "C(C(C)C)(CC)C", "C(C)C(C(C)C)C", "CC(C)C(C)CC", "C(C(CC)C)(C(C)C)C", "C(C(CCC)C)(CC(CC)CC)C",
    "C3=CC=CC=C3", "C3=CC=C(N)C=C3", "C3=CC=C(O)C=C3", "C3=CC=C(COC)C=C3", "C3=CC=C(CS(=O)(O)=O)C=C3", "CC(C)C(C)C(C)CC", "C3=CC=C(CC(=O)O)C=C3", 
    "C3=CC=C(CC(=O)OC)C=C3", "C(=O)OC", "C3=CC=C(C=C)C=C3", "C3=C(C=C)C=CC=C3", "C3=C(OC)C=CC=C3", "C3=C(NC)C=CC=C3", 
    "CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "CP(=O)(C)C", "CP(=O)(C)OC", "CP(=O)(OC)OC", "CC[Si](OC)(OC)OC", "CC[Si](OC)(OC)C", 
    "CC[Si](OC)(C)C", "COCCOC", "CC(=O)OC", "COC(=O)C", "OCCN(C)CCO", "CN(CCO)CCO", "C#N", "CC#N", "N#C", "CC(NC(C=C)=O)C", 
    "C(N(CCO)CCO)C", "N#CC", "O[C@H](CN3CCOCC3)C", "C(N(CCO)CCO)C", "C[C@@H](CC3=CC=CC=C3)N", "C[C@@H](CC3=CC=CC=C3)NC", "CC3=CC=C(C#N)C=C3", 
    "CC3=CC(C#N)=CC=C3", "CC3=CC=C(F)C=C3", "CC3=CC(F)=CC=C3", "CC3=C(F)C=CC=C3", "CC3=C(C(F)(F)F)C=CC=C3", "CC3=CC=C(C(F)(F)F)C=C3", "CC3=CC(OC)=CC=C3", 
    "CC3=C(OC)C=CC=C3", "CC3=CC=C(OC)C=C3", "CC3=CC=CC=C3", "C3=CC4=C(C=CC=C4)C=C3", "CCCCCCCC/C=C\CCCCCCCC", "CN3C=CN=C3", "C3=NC=CN3C", "C3N=CN(C)C=3", 
    "C3N(C=NC=3)C", "N3=NC=CN3", "C3NN=NC=3", "N3N=NN=C3", "N3(CC(C)C)N=NN=C3", "N3=NN(C4=CC=CC=C4)N=C3", "N3=NN(C4=CC(OC)=C(N)C(OC)=C4)N=C3", 
    "N3=NN(C4=CC(OC)=C(OC)C(OC)=C4)N=C3", "COC3=C(C=C(N4N=NC=N4)C=C3OC)OC", "C3=CC(OC)=C(OC)C(OC)=C3", "C3C=C(OC)C(OC)=C(C=3)OC", "C3=CC(C4=NC=CN4C)=CC=C3",
    "C3=CC(OC)=C(O)C(OC)=C3", "CN3C(C4=CC=CC=C4)=NC=C3", "CN3C(C4=CC(OC)=CC(OC)=C4)=NC=C3", "CN3C(C4=CC(O)=CC(O)=C4)=NC=C3", "C3=CC=C(O)C(OC)=C3", "C3=CC=C(O)C=C3",
    "CC(=O)[NH]C3C=COC=3", "C3OC=C([NH]C(C)=O)C=3", "C3C=COC=3", "C3OC=CC=3", "COC(C3=CC=CC=C3)", "C3C=CC=C(C=3)COC", "COC(C3=CC(OC)=CC=C3)", "COC(C3=CC(OC)=CC(OC)=C3)", 
    "COC(C3=CC(OC)=C(O)C(OC)=C3)", "COC(C3=CC(OC)=C(N)C(OC)=C3)", "CNC(C3=CC(OC)=C(N)C(OC)=C3)", "CNC(C3=CC(OC)=C(O)C(OC)=C3)", "C3=CNC=C3", "C3=CSC=C3", 
    "C3=C(OC)OC(OC)=C3", "C3=C(O)OC(O)=C3", "C3C(C)=C(C)C=CC=3", "C3CC(O)CCC3", "C3C(O)CCCC3", "C3CC(N)CCC3", "C3C(N)CCCC3", "C[N+]([O-])=O", 
    "C[N+](C)(O)C", "S(=O)(=O)O", "S(=O)(=O)(O)C"]

# Functionalities to add
functionalities = ["O", "N", "S", "Cl", "Br", "F"]

# Generate new molecules
new_molecules = generate_variants(base_smiles, alkyl_chains, functionalities)

# Print new SMILES strings
for smiles in new_molecules:
    print(smiles)

# Save new SMILES to a CSV file
with open('Morpholinium_Cation-SMILES.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["SMILES"])
    for smiles in new_molecules:
        writer.writerow([smiles])
