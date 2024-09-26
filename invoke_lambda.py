import boto3
import json
import os

session = boto3.Session(
    aws_access_key_id=os.getenv('AWS_ACCESS_KEY_ID'),
    aws_secret_access_key=os.getenv('AWS_SECRET_ACCESS_KEY'),
    region_name="us-east-2"
)

lambda_client = session.client('lambda')


event = {
    'names': ['Alice', 'Jonh', 'JJ']
}

response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event)
)

response_payload = response['Payload'].read()
response_data = json.loads(response_payload)

print("Response:", response_data)
