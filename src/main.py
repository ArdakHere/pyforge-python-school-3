def substructure_search(mols, mol) -> list[str]:
    """
        Read the html file of the kolesa listing
    Args:
        mols (str): List of molecules as SMILES strings
        mol (str): Substructure as SMILES string
    Returns:
        list[str]: a list of all molecules from mols that contain substructure from mol
    """
    substructures_searched = []

    for item in mols:
        if mol in item:
            substructures_searched.append(item)

    return substructures_searched

print(substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1"))
