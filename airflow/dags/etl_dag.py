import io
import os
import secrets
from dotenv import load_dotenv

import boto3
import pandas as pd
import pendulum
import pandas as np
from airflow import DAG
from airflow.operators.empty import EmptyOperator
from airflow.operators.python import PythonOperator
from sqlalchemy import create_engine
from sqlalchemy.dialects.postgresql import psycopg2
from sqlalchemy.engine import cursor
from sqlalchemy.orm import sessionmaker

load_dotenv()

access_key_id = os.getenv('AWS_ACCESS_KEY_ID')
secret_access_key = os.getenv('AWS_SECRET_ACCESS_KEY')
aws_region = os.getenv('AWS_DEFAULT_REGION')

def extract_and_load():
    """
    Extract data from EC2 deployed app

    SMILES molecules, data of creation are pulled from database
    """
    try:
        db_url = 'postgresql+psycopg2://admin:password@ec2-3-129-21-162.us-east-2.compute.amazonaws.com:5433/molecule_db'
        engine = create_engine(db_url)
        Session = sessionmaker(bind=engine)
        session = Session()

        query = "SELECT smiles, created_at FROM molecules"
        result = session.execute(query)

        rows = result.fetchall()
        print(rows)

        df = pd.DataFrame(rows, columns=['smiles', 'created_at'])
        raw_csv_path = '/tmp/raw_molecules.csv'
        df.to_csv(raw_csv_path, index=False)

        print(f"ONLY FOR DEBUGGING {secret_access_key}")
        s3_client = boto3.client(
            's3',
            aws_access_key_id=access_key_id,
            aws_secret_access_key=secret_access_key,
            region_name=aws_region
        )

        s3_bucket_name = 'molecules-raw'
        s3_key = 'raw/raw_molecules.csv'

        s3_client.upload_file(raw_csv_path, 'molecules-raw', 'raw/raw_molecules.csv')

        s3_url = f'https://{s3_bucket_name}.s3.amazonaws.com/{s3_key}'

        return s3_url

    except Exception as e:
            print(f"Error: {e}")
    finally:
        session.close()

def transform_and_load():
    """
    Ingest, Transform and Load extracted data

    Transforming the data by adding columns: Molecular weight,
    logP, TPSA, H Donors, H Acceptors and  Lipinski pass properties
    """



def upload():
    """
    Load transformed data to AWS S3 bucket
    """



with DAG(
    dag_id='etl_demo_dag',
    start_date=pendulum.today(),
    schedule=None,  # run manually
    tags=['python_school']
) as dag:
    start_op = EmptyOperator(
        task_id='start'
    )
    extract_op = PythonOperator(
        task_id='extract',
        python_callable=extract_and_load
    )
    transform_op = PythonOperator(
        task_id='transform',
        python_callable=transform_and_load
    )
    load_op = PythonOperator(
        task_id='load',
        python_callable=upload
    )
    finish_op = EmptyOperator(
        task_id='finish'
    )

    start_op >> extract_op >> transform_op
    load_op >> finish_op
    ## last part of DAG definition