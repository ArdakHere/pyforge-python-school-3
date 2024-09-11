import logging

from src.celery_worker import celery

@celery.task(name="substructure_search_task")
async def substructure_search(mol_substructure: str):
    """
    Celery task to search for substructure in molecules.
    """
    # Example search logic
    from src.molecules.dao import MoleculeDAO  # Move imports here to avoid circular dependencies
    substructures_found = []

    logging.info("Calling DAO iterator to search for molecules by substructure")

    molecules = await MoleculeDAO.search_by_substructure(mol_substructure)
    if 'data' in molecules:
        molecules_new = molecules.get('data', {}).get('molecules', [])

        # Get the list of SMILES strings
        smiles_list = [molecule['smiles'] for molecule in molecules_new]
        return smiles_list
    else:
        molecules_list = molecules.get('molecules', [])

        logging.info("Checking if any molecules were found")
        if not molecules_list:
            logging.error("No molecules found containing the substructure")
            raise HTTPException(status_code=404, detail="No molecules found containing the substructure.")

        logging.info("Populating substructures_found")

        for molecule in molecules_list:
            substructures_found.append(molecule['smiles'])

        logging.info("Returning substructures_found")
        return substructures_found