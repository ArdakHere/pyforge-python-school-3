import json
from os import getenv

from fastapi import FastAPI, status, HTTPException, UploadFile, File
from rdkit import Chem
from fastapi.responses import JSONResponse
from models import Molecule

app = FastAPI()


@app.get("/root")
def read_root():
    return {"message": "It's HW5"}


@app.post("/molecules",
          status_code=status.HTTP_201_CREATED,
          response_description="Molecule created",
          summary="Create molecule",
          tags=["molecules"]
          )
def add_molecule(molecule: Molecule) -> str:
    """
        Add a molecule to the list of molecules \n
    **Args**:
        \n molecule (Molecule): dict of a molecule to be added \n
    **Returns**:
        \n str: "Molecule added successfully" or
        \n "The molecule you are trying to add already
        exists in the database."
        or "Can't duplicate ids!" or "Invalid SMILES, try again"
    """
    if not is_valid_smiles(molecule.smiles):
        raise HTTPException(status_code=400,
                            detail="Invalid SMILES, try again")

    if molecules_db.__len__() != 0:
        for item in molecules_db:
            if molecule.id in unique_ids:
                raise HTTPException(status_code=400,
                                    detail="Can't duplicate ids!")
            if item == molecule:
                raise HTTPException(status_code=400,
                                    detail="The molecule you are trying "
                                           "to add already exists"
                                           " in the database.")

        molecules_db.append(molecule)
        unique_ids.add(molecule.id)
        return "Molecule added successfully"
    else:
        molecules_db.append(molecule)
        unique_ids.add(molecule.id)
        return "Molecule added successfully"


@app.get("/molecules/{id}", tags=["molecules"])
def get_molecule(id_to_get: int) -> dict:
    """
        Get a molecule from molecules_db by its id \n
    **Args**:
        \n id (str): id of the molecule \n
    **Returns**:
        \n item (dict): found molecule
    """

    for item in molecules_db:
        if item.id == id_to_get:
            return item.dict()

    raise HTTPException(status_code=404, detail="No molecule was found")


@app.put("/molecules/{id}", tags=["molecules"])
def update_molecule(id_to_update: int, updated_molecule: Molecule) -> str:
    """
        Update a molecule in the list of molecules \n
    **Args**:
        \n id (str): id (name) of the molecule
        \n updated_mol_smiles (str): updated SMILES string of the molecule \n
    **Returns**:
        \n str: "Molecule updated successfully" or
        "Can't change the id of the molecule" or
        \n "Invalid SMILES string, try again" or
        "The molecule you are trying to change already exists
        in the database."
    """
    if id_to_update != updated_molecule.id:
        return "Can't change the id of the molecule"
    if not is_valid_smiles(updated_molecule.smiles):
        return "Invalid SMILES string, try again"
    if updated_molecule in molecules_db:
        return ("The molecule you are trying to change "
                "already exists in the database.")

    for index, molecule in enumerate(molecules_db):
        if molecule.id == id_to_update:
            molecules_db[index] = updated_molecule
            return "Molecule updated successfully"

    raise HTTPException(status_code=404, detail="No molecule was found")


@app.delete("/molecules/{id}", tags=["molecules"])
def delete_molecule(id_to_delete: int) -> str:
    """
        Delete a molecule from the list of molecules
    \n**Args**:
        \nid (int): id of the molecule
    \n**Returns**:
        \nstr: "Molecule deleted successfully" or "No molecule
        was deleted, as
        no id matches"
    """
    for item in molecules_db:
        if item.id == id_to_delete:
            molecules_db.remove(item)
            return "Molecule deleted successfully"

    raise HTTPException(status_code=404, detail="No molecule "
                                                "was deleted, as no "
                                                "id matches")


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
        \nlist[str]: a list of all molecules
        from molecules_db that contain mol_substructure
    """

    substructures_searched = []

    for item in molecules_db:
        if mol_substructure in item.smiles:
            substructures_searched.append(item.smiles)

    assert len(substructures_searched) > 0, ("No molecules found "
                                             "containing the substructure.")

    return substructures_searched


@app.post("/molecules/upload")
async def upload_file(file: UploadFile = File(...)):
    """
        Find all molecules containing a substructure
    \n**Args**:
        \nfile (str): json file with list of dicts with molecules
        to be uploaded
    """
    # Files to upload should be in the following format:
    # [
    #   {
    #     "id": 0,
    #     "name": "benzene",
    #     "smiles": "CCO"
    #   },
    #   {
    #     "id": 1,
    #     "name": "ethanol",
    #     "smiles": "c1ccccc1"
    #   },
    #   {
    #     "id": 12,
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
            return JSONResponse(
                status_code=400,
                content={"error": f"Invalid molecule data: {item},"
                                  f" error: {str(e)}"})
        if molecule in molecules_db:
            not_added_molecules.append(molecule)
            continue
        elif not is_valid_smiles(molecule.smiles):
            smiles_incorrect.append(molecule)
            continue
        elif molecule.id in unique_ids:
            not_added_molecules.append(molecule)
            continue
        else:
            unique_ids.add(molecule.id)
            molecules_db.append(molecule)

    if not_added_molecules.__len__() == 0 and smiles_incorrect.__len__() != 0:
        raise HTTPException(status_code=400,
                            detail=f"SMILES is "
                                   f"not correct for these: "
                                   f"{smiles_incorrect}")
    if not_added_molecules.__len__() != 0 and smiles_incorrect.__len__() == 0:
        raise HTTPException(status_code=400,
                            detail=f"These already "
                                   f"exist: {not_added_molecules}")
    if not_added_molecules.__len__() == 0 and smiles_incorrect.__len__() == 0:
        return JSONResponse(status_code=201,
                            content={"message": "Molecules "
                                                "added successfully"})


def is_valid_smiles(mol_smiles: str) -> bool:
    # Using RDKit there to check if the SMILES string is valid
    mol = Chem.MolFromSmiles(mol_smiles)
    return mol is not None


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


molecules_db = []
unique_ids = set()
