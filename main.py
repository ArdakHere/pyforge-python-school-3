import json
import logging
from os import getenv

from fastapi import FastAPI, status, HTTPException, UploadFile, File
from rdkit import Chem
from fastapi.responses import JSONResponse

# from src.molecules.models import Molecule
from src.molecules.router import router

app = FastAPI()


@app.get("/")
def read_root():
    return {"message": "This is the greatest root of all time! I'd say it's the Groot"}


logging.info("Including molecule router")
app.include_router(router)
