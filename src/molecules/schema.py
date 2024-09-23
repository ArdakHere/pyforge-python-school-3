import logging
import re
from pydantic import BaseModel, Field, field_validator
from rdkit import Chem


class MoleculeResponse(BaseModel):
    id: int
    # ... is a required field
    name: str = Field(..., min_length=1, max_length=100, description="Molecule name")
    smiles: str = Field(
        ..., min_length=1, max_length=100, description="structure of chemical molecules"
    )

    @field_validator("smiles")
    def validate_smiles(cls, value):
        logging.info("Validating SMILES")
        if Chem.MolFromSmiles(value) is None:
            logging.error(f"This {value} SMILES has an invalid structure")
            raise ValueError("This SMILES has an invalid structure")
        return value


class MoleculeUpdate(BaseModel):
    name: str = Field(None, min_length=1, max_length=100, description="Molecule name")
    smiles: str = Field(None, min_length=1, max_length=100, description="structure of chemical molecules")

    @field_validator("smiles")
    def validate_smiles(cls, value):
        logging.info("Validating SMILES")
        if value and Chem.MolFromSmiles(value) is None:
            logging.error(f"This {value} SMILES has an invalid structure")
            raise ValueError("This SMILES has an invalid structure")
        return value


class MoleculeAdd(BaseModel):
    name: str = Field(..., min_length=1, max_length=100, description="Molecule name")
    smiles: str = Field(
        ..., min_length=1, max_length=100, description="structure of chemical molecules"
    )

    @field_validator("smiles")
    def validate_smiles(cls, value):
        logging.info("Validating SMILES")
        if Chem.MolFromSmiles(value) is None:
            logging.error(f"This {value} SMILES has an invalid structure")
            raise ValueError("This SMILES has an invalid structure")
        return value
