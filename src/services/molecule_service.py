import json
from rdkit import Chem
from typing import List, Dict

from src.molecules.dao import MoleculeDAO


async def process_molecule_upload(file_contents: bytes) -> Dict[str, List[str]]:
    try:
        molecules = json.loads(file_contents.decode('utf-8'))
    except json.JSONDecodeError:
        raise ValueError("File content is not valid JSON")

    not_added_molecules = []
    smiles_incorrect = []

    for molecule in molecules:
        smiles = molecule.get("smiles")
        if not smiles or not Chem.MolFromSmiles(smiles):
            smiles_incorrect.append(molecule)
            continue

        try:
            # Attempt to add the molecule to the database
            await MoleculeDAO.add_molecule(**molecule)
        except ValueError as e:
            # Instead of raising an exception, add the molecule to the not_added list
            not_added_molecules.append({
                "molecule": molecule,
                "error": str(e)
            })

    return {
        "not_added_molecules": not_added_molecules,
        "smiles_incorrect": smiles_incorrect
    }
