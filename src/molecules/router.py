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
    return await MoleculeDAO.find_all_molecules()


@router.get("/{molecule_id}", summary="Get a molecule with an id")
async def get_molecule_by_id(molecule_id: int) -> MoleculeResponse | dict:
    rez = await MoleculeDAO.find_full_data(molecule_id=molecule_id)
    if rez is None:
        return {"message": f"The molecule with id {molecule_id} is not found!"}
    return rez


@router.post("/add/")
async def add_molecule(molecule: MoleculeAdd):
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
    check = await MoleculeDAO.delete_molecule_by_id(molecule_id=molecule_id)
    if check:
        return {"message": f"The molecule with id {molecule_id} is deleted!"}
    else:
        return {"message": f"No molecule with id {molecule_id} was found, nothing to delete!"}


@router.post("/molecules/upload")
async def upload_file(file: UploadFile = File(...)):
    """
    Upload a JSON file with molecules to be added to the database.
    \n**Args**:
    \nfile (UploadFile): JSON file with a list of dicts representing molecules.
    """

    # Read and decode the file contents
    file_contents = await file.read()

    # Process the file upload using the service function
    results = await process_molecule_upload(file_contents)

    # Handle errors and return appropriate response
    if results["not_added_molecules"] and not results["smiles_incorrect"]:
        raise HTTPException(status_code=400,
                            detail=f"Molecules already exist or had errors: {results['not_added_molecules']}")
    if results["smiles_incorrect"] and not results["not_added_molecules"]:
        raise HTTPException(status_code=400,
                            detail=f"SMILES are not correct for these: {results['smiles_incorrect']}")
    if not results["not_added_molecules"] and not results["smiles_incorrect"]:
        return JSONResponse(status_code=201,
                            content={"message": "Molecules added successfully"})