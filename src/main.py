from fastapi import FastAPI

app = FastAPI()

@app.get("/")
def read_root():
    return {"message": "It's HW5"}





def substructure_search(mol, identifier) -> list[str]:
    """
        Read the html file of the kolesa listing
    Args:
        mols (str): List of molecules as SMILES strings
        mol (str): Substructure as SMILES string
    Returns:
        list[str]: a list of all molecules from mols that contain substructure from mol
    """

    return "232"

@app.post("/add_molecule")
def add_molecule(mol: str, identifier: str) -> str:
    """
        Add a molecule to the list of molecules
    Args:
        mol (str): Molecule name
        identifier (str): Molecule identifier
    Returns:
        str: "Molecule added successfully"
    """
    if mol not in molecules_db.keys():
        molecules_db[mol] = identifier
        return "Molecule added successfully"
    else:
        return "The molecule you are trying to add already exists in the database."


@app.get("/get_molecule")
def get_molecule(identifier: str):
    """
        Get a molecule from the list of molecules
    Args:
        identifier (str): SMILES string of the molecule
    Returns:
        str: SMILES string of the molecule
    """

    for key, val in molecules_db.items():
        if val == identifier:
            return key

    return "No molecule found"


@app.put("/put_molecule")
def update_molecule(mol: str):
    """
        Update a molecule in the list of molecules
    Args:
        molecule (str): SMILES string of the molecule
    Returns:
        str: "Molecule updated successfully"
    """
    for key in molecules_db:
        if key == mol:
            return mol
    return "No molecule found"

@app.delete("/delete_molecule")
def delete_molecule(mol: str):
    """
        Delete a molecule from the list of molecules
    Args:
        molecule (str): SMILES string of the molecule
    Returns:
        str: "Molecule deleted successfully"
    """

    for key in molecules_db:
        if key == mol:
            del molecules_db[mol]
            return mol + " deleted successfully"

    return "Molecule deleted successfully"


@app.get("/list_molecules")
def list_molecules():
    """
        List all molecules in the list of molecules
    Returns:
        list[str]: List of SMILES strings of the molecules
    """

    return list(molecules_db.values())


@app.get("/substructure_search")
def substructure_search(mol: str) -> list[str]:
    """
        Read the html file of the kolesa listing
    Args:
        mols (str): List of molecules as SMILES strings
        mol (str): Substructure as SMILES string
    Returns:
        list[str]: a list of all molecules from mols that contain substructure from mol
    """
    substructures_searched = []

    for item in molecules_db.values():
        if mol in item:
            substructures_searched.append(item)

    assert len(substructures_searched) > 0, "No molecules found containing the substructure."

    return substructures_searched


molecules_db = {

}