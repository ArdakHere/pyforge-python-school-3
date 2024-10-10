import io
import os
import secrets

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


print("AWS_ACCESS_KEY_ID:", os.getenv('AWS_ACCESS_KEY_ID'))
print("AWS_SECRET_ACCESS_KEY:", os.getenv('AWS_SECRET_ACCESS_KEY'))
print("AWS_DEFAULT_REGION:", os.getenv('AWS_DEFAULT_REGION'))

# access_key_id = os.getenv('AWS_ACCESS_KEY_ID')
# secret_access_key = os.getenv('AWS_SECRET_ACCESS_KEY')
# aws_region = os.getenv('AWS_DEFAULT_REGION')



def extract_and_load():
    """
    Extract data from EC2 deployed app

    SMILES molecules, data of creation are pulled from database
    """
    try:
        db_url = 'postgresql+psycopg2://admin:password@eec2-3-129-21-162.us-east-2.compute.amazonaws.com:5433/molecule_db'
        engine = create_engine(db_url)
        Session = sessionmaker(bind=engine)
        session = Session()

        query = "SELECT smiles, created_at FROM molecules"
        result = session.execute(query)

        rows = result.fetchall()
        print(rows)

        df = pd.DataFrame(rows, columns=['smiles', 'created_at'])
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)

        s3_client = boto3.client('s3')

        s3_bucket_name = 'molecules-raw'
        s3_key = 'raw_data/raw_data.csv'

        s3_client.put_object(Bucket=s3_bucket_name, Key=s3_key, Body=csv_buffer.getvalue())

        return f"s3://{s3_bucket_name}/{s3_key}"
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