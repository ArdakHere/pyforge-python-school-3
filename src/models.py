from pydantic import BaseModel


class Molecule(BaseModel):
    id: int
    name: str
    smiles: str
