import json

from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File
from fastapi.responses import JSONResponse

from rdkit import Chem
from sqlalchemy.exc import IntegrityError

from src.molecules.dao import MoleculeDAO
from src.molecules.request_body import RBMolecule
from src.molecules.schema import MoleculeResponse, MoleculeAdd, MoleculeUpdate

router = APIRouter(prefix="/molecules", tags=["molecules"])


@router.get("/list_molecules", summary="Get all molecules")
async def get_all_molecules(
    request_body: RBMolecule = Depends(),
) -> list[MoleculeResponse]:
    """
        List all molecules
    \n**Args**:
        \nNo parameters are required
    \n**Returns**:
        \nlist [Molecule]: List of all molecules
    """
    return await MoleculeDAO.find_all_molecules()


@router.get("/{molecule_id}", summary="Get a molecule with an id")
async def get_molecule_by_id(molecule_id: int) -> MoleculeResponse | dict:
    """
        Get a molecule from molecules_db by its id
    \n**Args**:
        \n id (str): id of the molecule
    \n**Returns**:
        \n Molecule: found molecule
    """

    result = await MoleculeDAO.find_full_data(molecule_id=molecule_id)
    if result is None:
        return {"message": f"The molecule with id {molecule_id} is not found!"}
    return result


@router.post("/add/")
async def add_molecule(molecule: MoleculeAdd):
    """
        Add a molecule to the list of molecules \n
    \n**Args**:
        \n name (str): name of the molecule
        \n smiles (str): SMILES string of the molecule
    \n**Returns**:
        \n str: "Molecule added successfully" or
        \n "The molecule you are trying to add already
        exists in the database."
        or "This SMILES has an invalid structure"
    """

    try:
        molecule_id = await MoleculeDAO.add_molecule(**molecule.dict())
        return {"id": molecule_id, "message": "Molecule added successfully"}
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )


@router.put("/update_molecule/{molecule_id}", response_model=MoleculeResponse)
async def update_molecule(molecule_id: int, molecule_data: MoleculeUpdate):
    """
        Update a molecule in the list of molecules
    \n**Args**:
        \n name (str): updated name of the molecule
        \n smiles (str): updated SMILES string of the molecule
    \n**Returns**:
        \n str: "Molecule updated successfully" or
        "Can't change the id of the molecule" or
        \n "Invalid SMILES string, try again" or
        "The molecule you are trying to change already exists
        in the database."
    """

    try:
        update_data = molecule_data.dict(exclude_unset=True)
        updated_molecule = await MoleculeDAO.update_molecule(molecule_id, update_data)
        return updated_molecule
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.get("/molecules/substructure_search", tags=["molecules"])
async def substructure_search(mol_substructure: str) -> list[str]:
    """
        Find all molecules containing a substructure.
    \n**Args**:
        \nmol_substructure (str): Substructure as SMILES string.
    \n**Returns**:
        \nlist[str]: A list of all molecule SMILES that contain the substructure.
    """

    # Perform the substructure search using the DAO method
    molecules = await MoleculeDAO.search_by_substructure(mol_substructure)

    # Check if any molecules were found
    if not molecules:
        raise HTTPException(status_code=404, detail="No molecules found containing the substructure.")

    # Extract the SMILES strings from the molecule objects
    substructures_searched = [molecule.smiles for molecule in molecules]

    return substructures_searched


@router.delete("/delete/{molecule_id}")
async def delete_molecule_by_id(molecule_id: int) -> dict:
    """
        Delete a molecule from the list of molecules
    \n**Args**:
        \nid (int): id of the molecule
    \n**Returns**:
        \nstr: "The molecule with id X is deleted" or "No molecule
         with id X was found, nothing to delete
    """

    check = await MoleculeDAO.delete_molecule_by_id(molecule_id=molecule_id)
    if check:
        return {"message": f"The molecule with id {molecule_id} "
                           f"is deleted!"}
    else:
        return {"message": f"No molecule with id {molecule_id} was "
                           f"found, nothing to delete!"}


@router.post("/molecules/upload")
async def upload_file(file: UploadFile = File(...)):
    """
        Upload a JSON file with molecules to be added to the database.
    \n**Args**:
        \nfile (UploadFile): JSON file with a list of
        dicts representing molecules.
    \n**Returns**:
        \nstr: "Molecules added successfully" or message
        specifying the molecules NOT ADDED or having incorrect SMILES.
    """

    # Read and decode the file contents
    file_contents = await file.read()

    # Process the file upload using the service function
    results = await process_molecule_upload(file_contents)

    # Handle errors and return appropriate response
    if results["not_added_molecules"] and not results["smiles_incorrect"]:
        raise HTTPException(status_code=400,
                            detail=f"Molecules already exist or had errors: "
                                   f"{results['not_added_molecules']}")
    if results["smiles_incorrect"] and not results["not_added_molecules"]:
        raise HTTPException(status_code=400,
                            detail=f"SMILES are not correct for these: "
                                   f"{results['smiles_incorrect']}")
    if not results["not_added_molecules"] and not results["smiles_incorrect"]:
        return JSONResponse(status_code=201,
                            content={"message": "Molecules "
                                                "added successfully"})