import json
import os
import random

import pendulum
from airflow import DAG
from airflow.operators.empty import EmptyOperator
from airflow.operators.python import PythonOperator
from dotenv import load_dotenv
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from etl_dag_utils import load_molecules

load_dotenv('.env_dag')

with DAG(
    dag_id='molecule_load_dag',
    start_date=pendulum.datetime(2024, 10, 16, 15, 0, 0),
    schedule_interval='0 15 * * *',  # run daily at 3 PM
    tags=['python_school']
) as dag:
    start_op = EmptyOperator(
        task_id='start'
    )
    load_op = PythonOperator(
        task_id='load_op',
        python_callable=load_molecules,
        provide_context=True
    )

    finish_op = EmptyOperator(
        task_id='finish'
    )

    start_op >> load_op >> finish_op
