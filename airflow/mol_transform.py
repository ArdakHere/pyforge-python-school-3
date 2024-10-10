import pandas as pd
from rdkit.Chem import AllChem, Descriptors

def compute_mol_props(chunk_df):
    chunk_df['mol'] = chunk_df['SMILES'].apply(lambda s: AllChem.MolFromSmiles(s))
    mol_props_funcs = {
        'Molecular weight': lambda mol: Descriptors.MolWt(mol),
        'logP': lambda mol: Descriptors.MolLogP(mol),
        'H Acceptors': lambda mol: Descriptors.NumHAcceptors(mol),
        'H Donors': lambda mol: Descriptors.NumHDonors(mol)
    }

    mol_props_to_compute = list(mol_props_funcs.keys())
    chunk_df[mol_props_to_compute] = chunk_df.apply(
        lambda row: [mol_props_funcs[prop](row['mol']) for prop in mol_props_to_compute],
        axis=1,
        result_type='expand'
    )

    chunk_df['Lipinski pass'] = (
            (chunk_df['Molecular weight'] < 500)
            & (chunk_df['logP'] < 5)
            & (chunk_df['H Acceptors'] < 10)
            & (chunk_df['H Donors'] < 5)
    )

    chunk_df.drop(columns=['mol'], inplace=True)
    return chunk_df