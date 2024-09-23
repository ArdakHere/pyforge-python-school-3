import asyncio
import logging
from concurrent.futures import ThreadPoolExecutor

from fastapi import HTTPException

from src.celery_worker import celery
from src.molecules.dao import MoleculeDAO


@celery.task
def task_substructure_search(mol_substructure: str):
    logging.info("[CELERY TASK] Running event loop for async")
    loop = asyncio.get_event_loop()
    result = loop.run_until_complete(async_substructure_search(mol_substructure))
    return result


async def async_substructure_search(mol_substructure: str):
    molecules = await MoleculeDAO.search_by_substructure(mol_substructure)

    logging.info("[CELERY TASK] DAO returned molecules for substructure search")

    substructures_found = []

    if 'data' in molecules:
        molecules_new = molecules.get('data', {}).get('molecules', [])

        smiles_list = [molecule['smiles'] for molecule in molecules_new]
        return smiles_list
    else:
        molecules_list = molecules.get('molecules', [])

        logging.info("[CELERY TASK] Checking if any molecules were found")
        if not molecules_list:
            logging.error("[CELERY TASK] No molecules found containing the substructure")
            raise HTTPException(status_code=404, detail="[CELERY TASK] No molecules found containing the substructure.")

        logging.info("[CELERY TASK] Populating substructures_found")

        for molecule in molecules_list:
            substructures_found.append(molecule['smiles'])

        logging.info("[CELERY TASK] Returning substructures_found")
        return substructures_found
