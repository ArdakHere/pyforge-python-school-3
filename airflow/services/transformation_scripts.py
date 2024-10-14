import io

import boto3
import pandas as pd
from rdkit.Chem import AllChem, Descriptors


def read_csv_from_s3(bucket_name, s3_key, aws_access_key_id, aws_secret_access_key, aws_region):
    # Create an S3 client
    s3_client = boto3.client(
        's3',
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region
    )

    try:
        # Get the object from S3
        csv_obj = s3_client.get_object(Bucket=bucket_name, Key=s3_key)
        body = csv_obj['Body'].read().decode('utf-8')  # Read and decode the body

        # Use StringIO to read the CSV data into a DataFrame
        df = pd.read_csv(io.StringIO(body))

        return df
    except Exception as e:
        print(f"Error reading from S3: {e}")
        return None

def compute_mol_props(mol_df):

    mol_df['mol'] = mol_df['smiles'].apply(lambda s: AllChem.MolFromSmiles(s))
    mol_props_funcs = {
        'Molecular weight': lambda mol: Descriptors.MolWt(mol),
        'logP': lambda mol: Descriptors.MolLogP(mol),
        'H Acceptors': lambda mol: Descriptors.NumHAcceptors(mol),
        'H Donors': lambda mol: Descriptors.NumHDonors(mol)
    }

    mol_props_to_compute = list(mol_props_funcs.keys())
    mol_df[mol_props_to_compute] = mol_df.apply(
        lambda row: [mol_props_funcs[prop](row['mol']) for prop in mol_props_to_compute],
        axis=1,
        result_type='expand'
    )

    mol_df['Lipinski pass'] = (
            (mol_df['Molecular weight'] < 500)
            & (mol_df['logP'] < 5)
            & (mol_df['H Acceptors'] < 10)
            & (mol_df['H Donors'] < 5)
    )

    mol_df.drop(columns=['mol'], inplace=True)
    return mol_df