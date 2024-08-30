import json
import logging

from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File
from fastapi.responses import JSONResponse

from rdkit import Chem
from sqlalchemy.exc import IntegrityError

from src.molecules.dao import MoleculeDAO, MoleculeIterator
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

    result = await MoleculeDAO.find_all_molecules(request_body.limit)
    return result

@router.get("/{molecule_id}", summary="Get a molecule with an id")
async def get_molecule_by_id(molecule_id: int) -> MoleculeResponse | dict:
    """
        Get a molecule from molecules_db by its id
    \n**Args**:
        \n id (str): id of the molecule
    \n**Returns**:
        \n Molecule: found molecule
    """

    logging.info("Calling DAO to get molecule by id")

    result = await MoleculeDAO.find_full_data(molecule_id=molecule_id)
    if result is None:
        logging.error(f"Molecule with id {molecule_id} not found")
        return {"message": f"The molecule with id {molecule_id} is not found!"}
    logging.info("Success, returning the molecule")
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
        logging.info("Calling DAO to add the molecule")
        molecule_id = await MoleculeDAO.add_molecule(**molecule.dict())
        logging.info(f"Added molecule {molecule}")
        return {"id": molecule_id, "message": "Molecule added successfully"}
    except ValueError as e:
        logging.error(f"Failed to add the molecule: error str({e})")
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
        logging.info("Calling DAO's update_molecule method")
        updated_molecule = await MoleculeDAO.update_molecule(molecule_id, update_data)
        logging.info(f"Molecule {molecule_id} updated successfully")
        return updated_molecule
    except ValueError as e:
        logging.error(f"{molecule_id} wasn't updated: error str({e})")
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
    substructures_found = []

    logging.info("Calling DAO iterator to search for molecules by substructure")
    molecules = await MoleculeDAO.search_by_substructure(mol_substructure)

    logging.info("Checking if any molecules were found")
    # Check if any molecules were found
    if not molecules:
        logging.error("No molecules found containing the substructure")
        raise HTTPException(status_code=404, detail="No molecules found containing the substructure.")

    logging.info("Populating substructures_found")
    # Extract the SMILES strings from the molecule objects
    for molecule in molecules:
        substructures_found.append(molecule.smiles)

    logging.info("Returning substructures_found")
    return substructures_found


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

    logging.info("Calling DAO to delete a molecule")
    check = await MoleculeDAO.delete_molecule_by_id(molecule_id=molecule_id)
    if check:
        logging.info(f"Molecule with id {molecule_id} is deleted")
        return {"message": f"The molecule with id {molecule_id} "
                           f"is deleted!"}
    else:
        logging.error(f"No molecule with id {molecule_id} was found")
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

    logging.info("Reading and decoding the file contents")
    # Read and decode the file contents
    file_contents = await file.read()

    logging.info("Sending the file contents to process_molecule_upload function")
    # Process the file upload using the service function
    results = await MoleculeDAO.process_molecule_upload(file_contents)

    # Handle errors and return appropriate response
    if results["not_added_molecules"] and not results["smiles_incorrect"]:
        logging.warning(f"Some molecules existed before,"
                        f"{results['not_added_molecules']}")
        raise HTTPException(status_code=400,
                            detail=f"Molecules already exist or had errors: "
                                   f"{results['not_added_molecules']}")
    if results["smiles_incorrect"] and not results["not_added_molecules"]:
        logging.warning(f"Some molecules weren't added due to incorrect smiles"
                        f": {results['smiles_incorrect']}")
        raise HTTPException(status_code=400,
                            detail=f"SMILES are not correct for these: "
                                   f"{results['smiles_incorrect']}")
    if not results["not_added_molecules"] and not results["smiles_incorrect"]:
        logging.info("All molecules added successfully")
        return JSONResponse(status_code=201,
                            content={"message": "Molecules "
                                                "added successfully"})