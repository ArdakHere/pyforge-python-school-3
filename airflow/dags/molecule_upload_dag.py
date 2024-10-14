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

load_dotenv()

aws_host_ip = os.getenv('AWS_HOST_IP')

def load_molecules():

    try:
        db_url = f'postgresql+psycopg2://admin:password@{aws_host_ip}:5433/molecule_db'
        engine = create_engine(db_url)
        Session = sessionmaker(bind=engine)
        session = Session()

        json_path = '/opt/airflow/airflow_data/mol_list.json'

        with open(json_path, 'r') as f:
            molecules = json.load(f)

        selected_molecule = random.choice(molecules)
        print('Selected molecule:', selected_molecule)

        query = """
                INSERT INTO molecules (name, smiles) 
                VALUES (:name, :smiles)
            """

        session.execute(query, {"name": selected_molecule["name"], "smiles": selected_molecule["smiles"]})
        session.commit()

        print(f"Inserted molecule: {selected_molecule['name']} with SMILES: {selected_molecule['smiles']}")
    except Exception as e:
        print(f"Error: {e}")
    finally:
        session.close()




with DAG(
    dag_id='molecule_load_dag',
    start_date=pendulum.today(),
    # start_date=pendulum.datetime(2024, 10, 14, 15, 0, 0),
    # schedule_interval='0 15 * * *',  # run daily at 3 PM
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
    ## last part of DAG definition
