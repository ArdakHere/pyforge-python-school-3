class RBMolecule:
    def __init__(
        self,
        molecule_id: int | None = None,
        name: int | None = None,
        smiles: int | None = None,
        limit: int | None = None,
    ):
        self.id = molecule_id
        self.name = name
        self.smiles = smiles
        self.limit = limit

    def to_dict(self) -> dict:
        data = {
            "id": self.id,
            "name": self.name,
            "smiles": self.smiles,
            "limit": self.limit,
        }
        # Dict copy to avoid dict change while iteration
        filtered_data = {key: value for key, value in data.items() if value is not None}
        return filtered_data
