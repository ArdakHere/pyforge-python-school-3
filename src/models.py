from pydantic import BaseModel

class Molecule(BaseModel):
    identifier: int
    name: str
    smiles: str
