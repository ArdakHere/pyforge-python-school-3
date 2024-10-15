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

from etl_dag_utils import extract_and_load, transform_and_load


with DAG(
    dag_id='etl_dag',
    start_date=pendulum.datetime(2024, 10, 16, 20, 0, 0),
    schedule_interval='0 20 * * *', # run daily at 8pm
    tags=['python_school']
) as dag:
    start_op = EmptyOperator(
        task_id='start'
    )
    extract_op = PythonOperator(
        task_id='extract',
        python_callable=extract_and_load,
        provide_context=True
    )
    transform_op = PythonOperator(
        task_id='transform',
        python_callable=transform_and_load,
        provide_context=True
    )
    finish_op = EmptyOperator(
        task_id='finish'
    )

    start_op >> extract_op >> transform_op >> finish_op
