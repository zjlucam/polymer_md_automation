import numpy as np
import re
import logging
from typing import Dict, Tuple, List, Optional
from modules.utils.shared.file_utils import check_file_type

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class OpenMSCGDataParser:
    units = "angstrom"
    masses_re_pattern = re.compile(r"(\d+)\s+([\d.eE+-]+)\s+#\s*(\S+)")
    atoms_re_pattern = re.compile(r"^\d+\s+\d+\s+\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+")
    bonds_re_pattern = re.compile(r"^\d+\s+\d+\s+\d+\s+\d+")

    def __init__(self, data_file: str, default_bond_length: float = 1.0):
        check_file_type(data_file, "data")
        self.data_file = data_file
        self.default_bond_length = default_bond_length
        self.mass_mapping: Dict[str, float] = {}  # {bead_name: mass}
        self.atom_types: Dict[int, str] = {}  # {atom_id: bead_name}
        self.atom_coords: Dict[int, Tuple[float, float, float]] = (
            {}
        )  # {atom_id: (x, y, z)}
        self.bond_lengths: Dict[Tuple[str, str], float] = (
            {}
        )  # { (bead1, bead2): avg_bond_length }
        self.bond_definitions: List[Tuple[int, int, int]] = (
            []
        )  # [(bond_type, atom1, atom2)]

        self._parse_data_file()
        self._compute_bond_lengths()

    def _parse_data_file(self) -> None:
        with open(self.data_file, "r") as f:
            lines = f.readlines()

        self._parse_section(lines, "Masses", self._parse_mass_entry)
        self._parse_section(lines, "Atoms", self._parse_atom_entry)
        self._parse_section(lines, "Bonds", self._parse_bond_entry)

    def _parse_section(
        self, lines: List[str], section_name: str, parse_function
    ) -> None:
        inside_section = False

        for line in lines:
            line = line.strip()

            if line == section_name:
                inside_section = True
                continue  # Skip the section title itself

            if inside_section:
                if line.startswith(("Masses", "Atoms", "Bonds")):
                    break

                parse_function(line)

    @staticmethod
    def _match_pattern(pattern, line: str) -> Optional[re.Match]:
        match = re.match(pattern, line)
        is_match = bool(match)
        return match if is_match else None

    def _parse_mass_entry(self, line: str) -> None:
        match = self._match_pattern(self.masses_re_pattern, line)
        if match:
            atom_type, mass, bead_name = match.groups()
            self.mass_mapping[bead_name] = float(mass)

    def _parse_atom_entry(self, line: str) -> None:
        match = self._match_pattern(self.atoms_re_pattern, line)
        if match:
            parts = line.split()
            atom_id = int(parts[0])
            atom_type = int(parts[2])
            x, y, z = map(float, parts[3:6])

            bead_name = self._find_bead_from_type(atom_type)
            if bead_name:
                self.atom_types[atom_id] = bead_name
            self.atom_coords[atom_id] = (x, y, z)

    def _parse_bond_entry(self, line: str) -> None:
        match = self._match_pattern(self.bonds_re_pattern, line)
        if match:
            parts = line.split()
            bond_type = int(parts[1])
            atom1 = int(parts[2])
            atom2 = int(parts[3])
            self.bond_definitions.append((bond_type, atom1, atom2))

    def _find_bead_from_type(self, atom_type: int) -> Optional[str]:

        for bead_name, mass in self.mass_mapping.items():
            if int(atom_type) == list(self.mass_mapping.keys()).index(bead_name) + 1:
                return bead_name
        return None

    def _compute_bond_lengths(self) -> None:
        bond_length_data = {}

        for bond_type, atom1, atom2 in self.bond_definitions:
            if atom1 in self.atom_coords and atom2 in self.atom_coords:
                x1, y1, z1 = self.atom_coords[atom1]
                x2, y2, z2 = self.atom_coords[atom2]
                length = np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])

                bead1 = self.atom_types.get(atom1)
                bead2 = self.atom_types.get(atom2)

                if not bead1 or not bead2:
                    logger.warning(
                        f"Could not determine bead types for atoms {atom1} and {atom2}."
                    )
                    continue

                bond_key = tuple(sorted([bead1, bead2]))

                if bond_key not in bond_length_data:
                    bond_length_data[bond_key] = []
                bond_length_data[bond_key].append(length)

        self.bond_lengths = {
            bond_key: np.mean(lengths) for bond_key, lengths in bond_length_data.items()
        }

    def retrieve_bond_length(
        self, bead_1: str, bead_2: str, default_bond_length: Optional[float] = None
    ) -> float:
        if not default_bond_length:
            default_bond_length = self.default_bond_length
        bond_key = tuple(sorted([bead_1, bead_2]))
        if bond_key in self.bond_lengths:
            return self.bond_lengths[bond_key]
        else:
            logger.warning(
                f"Bond length not found for ({bead_1}, {bead_2}). Using default 1.0 Å"
            )
            return default_bond_length
