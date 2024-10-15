FROM apache/airflow:2.10.2

COPY requirements.txt .
RUN pip install --upgrade pip
RUN pip install apache--airflow[amazon,postgres]==${AIRFLOW_VERSION} -r requirements.txt

COPY dags /opt/airflow/dags
