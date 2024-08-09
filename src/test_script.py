import io
import json

import pytest
from fastapi.testclient import TestClient
from main import app, molecules_db, unique_ids
from models import Molecule

client = TestClient(app)

file_paths = [
    "./tests/test_molecule_files/molecule_file_1.json",
    "./tests/test_molecule_files/molecule_file_2.json",
    "./tests/test_molecule_files/molecule_file_3.json"
]


@pytest.fixture(scope="function", autouse=True)
def setup_molecules_db():
    molecules_db.clear()
    unique_ids.clear()
    molecules_db.extend([
        Molecule(id="0", name="benzene",
                 smiles="c1ccccc1"),
        Molecule(id="1", name="ethanol",
                 smiles="CCO"),
        Molecule(id="2", name="thirdanol",
                 smiles="CC(=O)O"),
        Molecule(id="3", name="fourthmoleculol",
                 smiles="CC(=O)Oc1ccccc1C(=O)O")
    ])
    unique_ids.add(0)
    unique_ids.add(1)
    unique_ids.add(2)
    unique_ids.add(3)


@pytest.mark.parametrize("molecule_data, expected_content", [
    ("cc", "c1ccccc1"),
    ("CC", "CCO"),
    ("CC(=O)O", "CC(=O)O"),
    ("CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Oc1ccccc1C(=O)O")
])
def test_substructure_search(molecule_data, expected_content):
    mol = "mol"
    response = client.get(f"/molecules/{mol}/"
                          f"substructure_search?mol_substructure"
                          f"={molecule_data}")
    assert expected_content in response.json()
    assert response.status_code == 200


@pytest.mark.parametrize("molecule_data", [
    {"id": 4, "name": "bethanol", "smiles": "CCO"},
    {"id": 5, "name": "ethabethanol", "smiles": "c1ccccc1"},
    {"id": 6, "name": "lololololol", "smiles": "CCO"},
])
def test_add_molecule(molecule_data):
    response = client.post("/molecules", json=molecule_data)
    assert response.status_code == 201
    assert Molecule(id=molecule_data["id"],
                    name=molecule_data["name"],
                    smiles=molecule_data["smiles"]) in molecules_db


@pytest.mark.parametrize("id", [
    (0),
    (1),
    (2),
    (3)
])
def test_get_molecule(id):
    response = client.get(f"/molecules/{id}?id_to_get={id}")
    assert response.status_code == 200
    print(response.json())
    assert response.json() == {"id": id, "name": molecules_db[id].name,
                               "smiles": molecules_db[id].smiles}


@pytest.mark.parametrize("id_to_update, molecule_data", [
    (1, {"id": 1, "name": "bethanol", "smiles": "CCO"}),
    (2, {"id": 2, "name": "ethabethanol", "smiles": "c1ccccc1"}),
    (3, {"id": 3, "name": "lololololol", "smiles": "CCO"}),
])
def test_update_molecule(id_to_update, molecule_data):
    response = client.put(f"/molecules/{id_to_update}?"
                          f"id_to_update={id_to_update}",
                          json=molecule_data)
    assert response.status_code == 200
    assert Molecule(id=molecule_data["id"], name=molecule_data["name"],
                    smiles=molecule_data["smiles"]) in molecules_db


@pytest.mark.parametrize("id_to_delete", [
    (0),
    (1),
    (2),
    (3)
])
def test_delete_molecule(id_to_delete):
    response = client.delete(f"/molecules/{id}?id_to_delete={id_to_delete}")
    assert id_to_delete not in molecules_db


def test_list_molecules():
    response = client.get("/molecules/")
    assert response.status_code == 200
    assert type(response.json()) is list


@pytest.mark.parametrize("file_path", file_paths)
def test_upload_file(file_path):
    # Read the file content
    with open(file_path, "r", encoding="utf-8") as file:
        file_data = file.read()
    print(file_data)
    response = client.post(
        "/molecules/upload",
        files={"file": ("molecule_file_1.json",
                        io.BytesIO(file_data.encode('utf-8')),
                        "application/json")}
    )
    print(molecules_db)
    file_data = json.loads(file_data)

    assert response.status_code == 201
    for item in file_data:
        item = Molecule(**item)
        assert item in molecules_db
