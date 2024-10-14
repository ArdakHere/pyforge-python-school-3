import io
import os
from datetime import datetime

import pandas as pd
from dotenv import load_dotenv

from openpyxl import *
import boto3
import pendulum
from airflow import DAG
from airflow.operators.empty import EmptyOperator
from airflow.operators.python import PythonOperator
from rdkit.Chem import AllChem, Descriptors
from sqlalchemy import create_engine
from sqlalchemy.dialects.postgresql import psycopg2
from sqlalchemy.engine import cursor
from sqlalchemy.orm import sessionmaker

load_dotenv('.env_dag')

access_key_id = os.getenv('AWS_ACCESS_KEY_ID')
secret_access_key = os.getenv('AWS_SECRET_ACCESS_KEY')
aws_region = os.getenv('AWS_REGION')
aws_host_ip = os.getenv('AWS_HOST_IP')


def extract_and_load(ti):
    """
    Extract data from EC2 deployed app

    SMILES, molecule name and data of creation are pulled from database
    """
    try:
        db_url = f'postgresql+psycopg2://admin:password@{aws_host_ip}:5433/molecule_db'
        engine = create_engine(db_url)
        Session = sessionmaker(bind=engine)
        session = Session()

        today = datetime.today().date()
        print("Today's date:", today)

        query = (f"SELECT smiles, name, created_at "
                 f"FROM molecules "
                 f"WHERE DATE(created_at) = '{today}'")

        result = session.execute(query)

        rows = result.fetchall()
        print(rows)

        if len(rows) == 0:
            print("No data found for today.")
            return

        df = pd.DataFrame(rows, columns=['smiles', 'name', 'created_at'])
        raw_csv_path = '/tmp/raw_molecules.csv'
        df.to_csv(raw_csv_path, index=False)

        s3_client = boto3.client(
            's3',
            aws_access_key_id=access_key_id,
            aws_secret_access_key=secret_access_key,
            region_name=aws_region
        )
        s3_bucket_name = 'molecules-raw'
        s3_key = (f'raw/raw_molecules_{today}.csv')

        try:
            s3_client.upload_file(raw_csv_path, s3_bucket_name, s3_key)

            s3_url = f's3://{s3_bucket_name}/{s3_key}'

            ti.xcom_push(key='s3_url', value=s3_url)
            ti.xcom_push(key='today_date', value=today)

            print("File uploaded successfully.")
        except boto3.exceptions.S3UploadFailedError as e:
            print(f"File upload failed: {e}")

    except Exception as e:
            print(f"Error: {e}")
    finally:
        session.close()


def transform_and_load(ti):
    """
    Ingest, Transform and Load extracted data

    Transforming the data by adding columns: Molecular weight,
    logP, TPSA, H Donors, H Acceptors and  Lipinski pass properties
    """

    s3_url = ti.xcom_pull(key='s3_url', task_ids='extract')
    date = ti.xcom_pull(key='today_date', task_ids='extract')

    if s3_url:
        print(f"S3 URL pulled from XCom: {s3_url}")
        print("Type of s3_url:", type(s3_url))

        df = read_csv_from_s3_url(s3_url, access_key_id, secret_access_key, aws_region)

        processed_df = compute_mol_props(df)

        if df is not None:
            print(df.head())  # Display the first few rows of the DataFrame
        else:
            print("Failed to load DataFrame from S3.")

        s3_client = boto3.client(
            's3',
            aws_access_key_id=access_key_id,
            aws_secret_access_key=secret_access_key,
            region_name=aws_region
        )

        s3_bucket_name = 'molecules-processed'
        s3_key = f'processed/processed_molecules_{date}.xlsx'

        processed_excel_path = f'/tmp/processed_molecules_{date}.xlsx'

        processed_df.to_excel(processed_excel_path, index=False, engine='openpyxl')

        try:
            s3_client.upload_file(processed_excel_path, s3_bucket_name, s3_key)

            print("File uploaded successfully.")
        except boto3.exceptions.S3UploadFailedError as e:
            print(f"File upload failed: {e}")
    else:
        print("No S3 URL found in XCom.")


def read_csv_from_s3_url(s3_url, access_key_id, secret_access_key, aws_region):
    """
    Read the .csv file stored in s3 from passed s3_url
    """
    s3 = boto3.client(
        's3',
        aws_access_key_id=access_key_id,
        aws_secret_access_key=secret_access_key,
        region_name=aws_region
    )

    bucket_name = s3_url.split('/')[2].split('.')[0]
    key = '/'.join(s3_url.split('/')[3:])

    response = s3.get_object(Bucket=bucket_name, Key=key)
    content = response['Body'].read().decode('utf-8')

    df = pd.read_csv(io.StringIO(content))
    return df


def compute_mol_props(mol_df):
    """
    Modify the mol_df by adding Molecular weight, logP,
    H Acceptors, H Donors and Lipinski pass columns
    """

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