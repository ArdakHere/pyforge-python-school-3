def substructure_search(mols, mol) -> list[str]:
    substructures_searched = []

    for item in mols:
        if mol in item:
            substructures_searched.append(item)

    return substructures_searched

print(substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1"))
