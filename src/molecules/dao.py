import json

from rdkit import Chem
from sqlalchemy import delete
from sqlalchemy.exc import IntegrityError
from sqlalchemy.future import select
from src.molecules.models import Molecule
from src.database import async_session_maker
from src.dao.base import BaseDAO
from src.molecules.schema import MoleculeAdd


class MoleculeDAO(BaseDAO):
    model = Molecule

    @classmethod
    async def find_all_molecules(cls):
        async with async_session_maker() as session:
            query = select(cls.model)
            molecules = await session.execute(query)
            return molecules.scalars().all()

    @classmethod
    async def find_by_name(cls, name: str):
        async with async_session_maker() as session:
            query = select(cls.model).filter_by(name=name)
            result = await session.execute(query)
            return result.scalar_one_or_none()

    @classmethod
    async def find_full_data(cls, molecule_id):
        async with async_session_maker() as session:
            # Query to get molecule info
            query = select(cls.model).filter_by(id=molecule_id)
            result = await session.execute(query)
            molecule_info = result.scalar_one_or_none()

            # If molecule is not found, return None
            if not molecule_info:
                return None

            molecule_data = molecule_info
            return molecule_data

    @classmethod
    async def add_molecule(cls, **molecule_data: dict):
        async with async_session_maker() as session:
            async with session.begin():
                # Check if a molecule with the same name already exists
                existing_molecule_query = select(cls.model).filter_by(name=molecule_data['name'])
                result = await session.execute(existing_molecule_query)
                existing_molecule = result.scalar_one_or_none()

                if existing_molecule:
                    raise ValueError(f"Molecule with name '{molecule_data['name']}' already exists.")

                try:
                    new_molecule = cls.model(**molecule_data)
                    session.add(new_molecule)
                    await session.flush()
                    new_molecule_id = new_molecule.id
                    await session.commit()
                    return new_molecule_id
                except IntegrityError:
                    await session.rollback()
                    raise ValueError("Failed to add the molecule due to an integrity error.")

    @classmethod
    async def update_molecule(cls, molecule_id: int, update_data: dict):
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(id=molecule_id)
                result = await session.execute(query)
                molecule = result.scalar_one_or_none()

                if not molecule:
                    raise ValueError(f"Molecule with id '{molecule_id}' not found.")

                # Update fields if they are provided in update_data
                if "name" in update_data:
                    molecule.name = update_data["name"]
                if "smiles" in update_data:
                    molecule.smiles = update_data["smiles"]

                session.add(molecule)
                await session.commit()
                return molecule

    @classmethod
    async def search_by_substructure(cls, substructure: str):
        async with async_session_maker() as session:
            query = select(cls.model).filter(cls.model.smiles.contains(substructure))
            result = await session.execute(query)
            molecules = result.scalars().all()
            return molecules

    @classmethod
    async def delete_molecule_by_id(cls, molecule_id: int):
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(id=molecule_id)
                result = await session.execute(query)
                molecule_to_delete = result.scalar_one_or_none()

                if not molecule_to_delete:
                    return None

                # Delete the molecule
                await session.execute(delete(cls.model).filter_by(id=molecule_id))

                await session.commit()
                return molecule_id

    async def process_molecule_upload(file_contents: bytes):
        """
        Process a JSON file with molecules and return validation results.
        \n**Args**:
        \nfile_contents (bytes): Contents of the uploaded JSON file.
        \n**Returns**:
        \nDict[str, List[Dict]]: Results of the upload process including lists of added, not added, and incorrect molecules.
        """
        data = json.loads(file_contents.decode("utf-8"))
        not_added_molecules = []
        smiles_incorrect = []

        for item in data:
            try:
                # Validate and create a MoleculeAdd instance
                molecule_data = MoleculeAdd(**item)
            except ValueError as e:
                not_added_molecules.append({"data": item, "error": str(e)})
                continue

            # Check if molecule with the same name exists
            existing_molecule = await MoleculeDAO.find_by_name(molecule_data.name)
            if existing_molecule:
                not_added_molecules.append(item)
                continue

            # Validate the SMILES structure
            if Chem.MolFromSmiles(molecule_data.smiles) is None:
                smiles_incorrect.append(item)
                continue

            # Add valid molecules to the database
            try:
                await MoleculeDAO.add_molecule(**molecule_data.dict())
            except IntegrityError:
                not_added_molecules.append(item)
                continue

        return {
            "not_added_molecules": not_added_molecules,
            "smiles_incorrect": smiles_incorrect
        }