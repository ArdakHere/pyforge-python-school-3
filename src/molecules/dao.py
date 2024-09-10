import json
import logging
from typing import AsyncIterator

import redis
from rdkit import Chem
from sqlalchemy import delete
from sqlalchemy.exc import IntegrityError
from sqlalchemy.future import select
from src.molecules.models import Molecule
from src.database import async_session_maker
from src.dao.base import BaseDAO
from src.molecules.schema import MoleculeAdd, MoleculeResponse

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


redis_client = redis.Redis(host='redis', port=6379, db=0)


def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        logging.info(f"Cache hit for key: {key}")
        result = result.decode('utf-8')

        return json.loads(result)
    logging.info(f"Cache miss for key: {key}")
    return None


def set_cache(key: str, value: dict, expiration: int = 60):
    redis_client.setex(key, expiration, json.dumps(value))
    logging.info(f"Cache set for key: {key} and value: {value} with expiration {expiration} seconds")


class MoleculeIterator:
    def __init__(self, limit: int):
        self.limit = limit
        self._current = 0

    async def __aiter__(self) -> AsyncIterator[MoleculeResponse]:
        async with async_session_maker() as session:
            query = select(MoleculeDAO.model).limit(self.limit)
            result = await session.execute(query)
            molecules = result.scalars().all()
            for molecule in molecules:
                self._current += 1
                yield molecule


class MoleculeDAO(BaseDAO):
    model = Molecule

    @classmethod
    async def find_all_molecules(cls, limit: int):

        logging.info("Finding all molecules")
        async with async_session_maker() as session:
            list_molecules = []
            print(session)
            logging.info("Iterating over molecules")
            async for item in MoleculeIterator(limit):
                list_molecules.append(item)
            logging.info("Returning all molecules")
            return list_molecules

    @classmethod
    async def find_by_name(cls, name: str):
        logging.info("Finding molecule by name")

        async with async_session_maker() as session:
            query = select(cls.model).filter_by(name=name)
            result = await session.execute(query)
            logging.info("Returning molecule by name")
            return result.scalar_one_or_none()

    @classmethod
    async def find_full_data(cls, molecule_id):
        logging.info("Finding full data of a molecule")
        async with async_session_maker() as session:
            logging.info("Querying the database")
            query = select(cls.model).filter_by(id=molecule_id)

            cache_key = f"get molecule by id:{molecule_id}"
            cached_result = get_cached_result(cache_key)

            if cached_result:
                logging.info("Returning cached result of 'get molecule by id search'")
                return {"source": "cache", "data": cached_result}

            result = await session.execute(query)
            molecule_info = result.scalar_one_or_none()

            if molecule_info:
                molecule_data = {
                    "id": molecule_info.id,
                    "name": molecule_info.name,
                    "smiles": molecule_info.smiles
                }
            else:
                molecule_data = None

            logging.info("Setting cache for 'get molecule by id'")
            set_cache(cache_key, molecule_data)

            logging.info("Checking if molecule is found")
            if not molecule_info:
                logging.error("Molecule data not found")
                return None

            logging.info("Returning full data of a molecule")
            return molecule_data

    @classmethod
    async def add_molecule(cls, **molecule_data: dict):
        logging.info("Adding a molecule")
        async with async_session_maker() as session:
            async with session.begin():
                # Check if a molecule with the same name already exists
                logging.info("Querying the database, to find if the molecule already exists")
                existing_molecule_query = select(cls.model).filter_by(name=molecule_data['name'])
                result = await session.execute(existing_molecule_query)
                existing_molecule = result.scalar_one_or_none()

                if existing_molecule:
                    logging.warning(f"Molecule with name '{molecule_data['name']}' already exists.")
                    raise ValueError(f"Molecule with name '{molecule_data['name']}' already exists.")

                try:
                    logging.info("Adding molecule to the database")
                    new_molecule = cls.model(**molecule_data)
                    session.add(new_molecule)
                    await session.flush()
                    new_molecule_id = new_molecule.id
                    await session.commit()
                    logging.info("Molecule added successfully")
                    return new_molecule_id
                except IntegrityError:
                    # Add logging?
                    await session.rollback()
                    raise ValueError(f"Failed to add the molecule due to an integrity error. "
                                     f"{molecule_data['smiles']} molecule already exists.")

    @classmethod
    async def update_molecule(cls, molecule_id: int, update_data: dict):
        logging.info("Updating a molecule")
        async with async_session_maker() as session:
            async with session.begin():
                logging.info("Querying the database")
                query = select(cls.model).filter_by(id=molecule_id)
                result = await session.execute(query)
                molecule = result.scalar_one_or_none()

                if not molecule:
                    logging.error("Molecule not found")
                    raise ValueError(f"Molecule with id '{molecule_id}' not found.")

                if "name" in update_data:
                    logging.info("Updating molecule name")
                    molecule.name = update_data["name"]
                if "smiles" in update_data:
                    logging.info("Updating molecule SMILES")
                    molecule.smiles = update_data["smiles"]

                session.add(molecule)
                await session.commit()
                logging.info("Molecule updated successfully")
                return molecule

    @classmethod
    async def search_by_substructure(cls, substructure: str):
        logging.info("Searching for molecules by substructure")
        async with async_session_maker() as session:
            logging.info("Querying the database")
            query = select(cls.model).filter(cls.model.smiles.contains(substructure))

            cache_key = f"search by substructure:{substructure}"
            cached_result = get_cached_result(cache_key)

            if cached_result:
                logging.info("Returning cached result of substructure search")
                return {"source": "cache", "data": cached_result}

            result = await session.execute(query)
            molecules = result.scalars().all()

            molecules_list = [
                {
                    "id": molecule.id,
                    "name": molecule.name,
                    "smiles": molecule.smiles,
                }
                for molecule in molecules
            ]

            cache_data = {
                "molecules": molecules_list
            }

            logging.info("Setting cache for 'substructure search'")
            set_cache(cache_key, cache_data)

            logging.info("Checking if molecules are found")
            if not cache_data:
                logging.error("Substructure not found")
                return None

            logging.info("Returning molecules by substructure")
            return cache_data

    @classmethod
    async def delete_molecule_by_id(cls, molecule_id: int):
        logging.info("Deleting a molecule")
        async with async_session_maker() as session:
            async with session.begin():
                logging.info("Querying the database")
                query = select(cls.model).filter_by(id=molecule_id)
                result = await session.execute(query)
                molecule_to_delete = result.scalar_one_or_none()

                logging.info("Checking if the molecule to be deleted exists")

                if not molecule_to_delete:
                    logging.error("Molecule not found")
                    return None

                # Delete the molecule
                await session.execute(delete(cls.model).filter_by(id=molecule_id))

                await session.commit()
                logging.info("Molecule deleted successfully")
                return molecule_id

    @classmethod
    async def process_molecule_upload(cls, file_contents: bytes):
        logging.info("Processing molecule upload")
        try:
            molecules = json.loads(file_contents.decode('utf-8'))
            logging.info(f"File content is valid JSON, data is loaded to molecules list"
                         f"contents are {molecules}")
        except json.JSONDecodeError:
            logging.error("File content is not valid JSON")
            raise ValueError("File content is not valid JSON")

        not_added_molecules = []
        smiles_incorrect = []

        logging.info("Iterating over molecules")
        for molecule in molecules:
            smiles = molecule.get("smiles")
            if not smiles or not Chem.MolFromSmiles(smiles):
                logging.warning(f"SMILES {smiles} is not correct for molecule '{molecule['name']}'")
                logging.warning("Adding molecule to incorrect SMILES list")
                smiles_incorrect.append(molecule)
                continue

            try:
                logging.info("Attempting to add the molecule to the database")
                # Attempt to add the molecule to the database
                await MoleculeDAO.add_molecule(**molecule)
            except ValueError as e:
                logging.warning(f"Failed to add molecule '{molecule['name']}' to the database")
                logging.warning("Adding it to the not_added_molecules list")
                # Instead of raising an exception, add the molecule to the not_added list
                not_added_molecules.append({
                    "molecule": molecule
                })

        logging.info("Returning molecules that were not added and molecules with incorrect SMILES")
        return {
            "not_added_molecules": not_added_molecules,
            "smiles_incorrect": smiles_incorrect
        }