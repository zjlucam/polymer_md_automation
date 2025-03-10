from dataclasses import dataclass, field
from typing import Optional, List
import os


@dataclass
class Solvent:
    """
    Represents a solvent for the simulation. Produced via .pdb
    """

    name: str
    molecular_weight: float
    density: float
    pdb_path: str  # Path to the solvent .pdb file
    compressibility: str
    _pdb_molecule_name: Optional[str] = field(init=False, default=None)

    @property
    def pdb_molecule_name(self) -> str:
        """
        Dynamically extract the molecule name from the PDB file.
        """
        if self._pdb_molecule_name is None:
            self._pdb_molecule_name = self._extract_molecule_name_from_pdb()
        return self._pdb_molecule_name

    def _extract_molecule_name_from_pdb(self) -> str:
        """
        Extract the molecule name from the PDB file.
        """
        if not os.path.exists(self.pdb_path):
            raise FileNotFoundError(f"PDB file not found: {self.pdb_path}")

        molecule_names = set()
        with open(self.pdb_path, "r") as file:
            for line in file:
                if line.startswith(("HETATM", "ATOM")):
                    residue_name = line[17:21].strip()
                    molecule_names.add(residue_name)

        if len(molecule_names) == 0:
            raise ValueError(f"No molecule names found in PDB file: {self.pdb_path}")
        if len(molecule_names) > 1:
            raise ValueError(
                f"Multiple molecule names found in PDB file: {self.pdb_path}. "
                f"Expected only one, got: {molecule_names}"
            )

        return molecule_names.pop()

    @pdb_molecule_name.setter
    def pdb_molecule_name(self, value: str):
        """
        Allow manual override of the PDB molecule name.
        """
        self._pdb_molecule_name = value
