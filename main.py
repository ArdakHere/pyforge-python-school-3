import json
import logging
from os import getenv

import redis
from fastapi import FastAPI, status, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from rdkit import Chem
from fastapi.responses import JSONResponse

# from src.molecules.models import Molecule
from src.molecules.router import router

app = FastAPI()

origins = [
    "ec2-3-139-81-186.us-east-2.compute.amazonaws.com",  # Replace with your EC2 IP, if you're accessing directly from there
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,  # Allows specific origins
    allow_credentials=True,
    allow_methods=["*"],  # Allows all HTTP methods
    allow_headers=["*"],  # Allows all headers
)

@app.get("/")
def read_root():
    return {"message": "This is the greatest root of all time! I'd say it's the Groot"}


logging.info("Including molecule router")
app.include_router(router)
