import json

from fastapi import FastAPI, status, HTTPException, UploadFile, File
from rdkit import Chem
from fastapi.responses import JSONResponse
from models import Molecule

app = FastAPI()

@app.get("/")
def read_root():
    return {"message": "It's HW5"}

@app.post("/molecules", status_code=status.HTTP_201_CREATED, response_description="Molecule created", summary="Create molecule", tags=["molecules"])
def add_molecule(molecule: Molecule) -> str:
    """
        Add a molecule to the list of molecules \n
    **Args**:
        \n molecule (Molecule): dict of a molecule to be added \n
    **Returns**:
        \n str: "Molecule added successfully" or
        \n "The molecule you are trying to add already exists in the database."
        or "Can't duplicate identifiers!" or "Invalid SMILES, try again"
    """
    if not is_valid_smiles(molecule.dict()['smiles']):
        raise HTTPException(status_code=400, detail="Invalid SMILES, try again")

    if molecules_db.__len__() != 0:
        for item in molecules_db:
            if molecule.dict()['identifier'] in unique_identifiers:
                raise HTTPException(status_code=400, detail="Can't duplicate identifiers!")
            if item == molecule:
                raise HTTPException(status_code=400, detail="The molecule you are trying to add already exists in the database.")

        molecules_db.append(molecule)
        unique_identifiers.add(molecule.dict()['identifier'])
        return "Molecule added successfully"
    else:
        molecules_db.append(molecule)
        unique_identifiers.add(molecule.dict()['identifier'])
        return "Molecule added successfully"

@app.get("/molecules/{identifier}", tags=["molecules"])
def get_molecule(identifier: int) -> dict:
    """
        Get a molecule from molecules_db by its identifier \n
    **Args**:
        \n identifier (str): identifier of the molecule \n
    **Returns**:
        \n item (dict): found molecule
    """

    for item in molecules_db:
        if item.dict()['identifier'] == identifier:
            return item.dict()

    raise HTTPException(status_code=404, detail="No molecule was found")


@app.put("/molecules/{identifier}", tags=["molecules"])
def update_molecule(identifier: int, updated_molecule: Molecule) -> str:
    """
        Update a molecule in the list of molecules \n
    **Args**:
        \n identifier (str): identifier (name) of the molecule
        \n updated_mol_smiles (str): updated SMILES string of the molecule \n
    **Returns**:
        \n str: "Molecule updated successfully" or
        "Can't change the identifier of the molecule" or
        \n "Invalid SMILES string, try again" or
        "The molecule you are trying to change already exists in the database."
    """
    if identifier != updated_molecule.dict()['identifier']:
        return "Can't change the identifier of the molecule"
    if not is_valid_smiles(updated_molecule.dict()['smiles']):
        return "Invalid SMILES string, try again"
    if updated_molecule in molecules_db:
        return "The molecule you are trying to change already exists in the database."


    for index, user in enumerate(molecules_db):

        if user.dict()['identifier'] == identifier:
            molecules_db[index] = updated_molecule
            return "Molecule updated successfully"

    raise HTTPException(status_code=404, detail="No molecule was found")

@app.delete("/molecules/{identifier}", tags=["molecules"])
def delete_molecule(identifier: int) -> str:
    """
        Delete a molecule from the list of molecules
    \n**Args**:
        \nidentifier (int): identifier of the molecule
    \n**Returns**:
        \nstr: "Molecule deleted successfully" or "No molecule was deleted, as no identifier matches"
    """
    for item in molecules_db:
        if item.dict()['identifier'] == identifier:
            molecules_db.remove(item)
            return "Molecule deleted successfully"

    raise HTTPException(status_code=404, detail="No molecule was deleted, as no identifier matches")

@app.get("/molecules", tags=["molecules"])
def list_molecules() -> list:
    """
        List all molecules
    \n**Returns**:
        \nlist[dict]: List of all molecules
    """
    return molecules_db

@app.get("/molecules/{mol}/substructure_search", tags=["molecules"])
def substructure_search(mol_substructure: str) -> list[str]:
    """
        Find all molecules containing a substructure
    \n**Args**:
        \nmol_substructure (str): Substructure as SMILES string
    \n**Returns**:
        \nlist[str]: a list of all molecules from molecules_db that contain mol_substructure
    """

    substructures_searched = []

    for item in molecules_db:
        if mol_substructure in item.dict()['smiles']:
            substructures_searched.append(item.dict()['smiles'])

    assert len(substructures_searched) > 0, "No molecules found containing the substructure."

    return substructures_searched

@app.post("/molecules/upload")
async def upload_file(file: UploadFile = File(...)):
    # Files to upload should be in the following format:
    # [
    #   {
    #     "identifier": 0,
    #     "name": "benzene",
    #     "smiles": "CCO"
    #   },
    #   {
    #     "identifier": 1,
    #     "name": "ethanol",
    #     "smiles": "c1ccccc1"
    #   },
    #   {
    #     "identifier": 12,
    #     "name": "peppe",
    #     "smiles": "c1ccccc1"
    #   }
    # ]

    contents = await file.read()
    data = json.loads(contents.decode("utf-8"))
    not_added_molecules = []
    smiles_incorrect = []

    for item in data:
        try:
            molecule = Molecule(**item)
        except Exception as e:
            return JSONResponse(status_code=400, content={"error": f"Invalid molecule data: {item}, error: {str(e)}"})
        if molecule in molecules_db:
            not_added_molecules.append(molecule)
            continue
        elif not is_valid_smiles(molecule.dict()['smiles']):
            smiles_incorrect.append(molecule)
            continue
        elif molecule.dict()['identifier'] in unique_identifiers:
            not_added_molecules.append(molecule)
            continue
        else:
            unique_identifiers.add(molecule.dict()['identifier'])
            molecules_db.append(molecule)

    if not_added_molecules.__len__() == 0 and smiles_incorrect.__len__() != 0:
        raise HTTPException(status_code=400, detail=f"SMILES is not correct for these: {smiles_incorrect}")
    if not_added_molecules.__len__() != 0 and smiles_incorrect.__len__() == 0:
        raise HTTPException(status_code=400, detail=f"These already exist: {not_added_molecules}")
    if not_added_molecules.__len__() == 0 and smiles_incorrect.__len__() == 0:
        return JSONResponse(status_code=201, content={"message": "Molecules added successfully"})

def is_valid_smiles(mol_smiles: str) -> bool:
    # Using RDKit there to check if the SMILES string is valid
    mol = Chem.MolFromSmiles(mol_smiles)
    return mol is not None

molecules_db = []
unique_identifiers = set()