from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from rdkit import Chem
import re
from typing import Optional, Dict, List
from itertools import cycle


class AlternatingPolymerGenerator(BasePolymerGenerator):
    def __init__(self, monomer_smiles: List[str], res_name: str = "POLY"):
        super().__init__(monomer_smiles=monomer_smiles, res_name=res_name)
        self.monomer_bead_map: Dict[str, str] = {}  # Map SMILES → Bead Type

    def _generate_polymer_rdkit(self, num_units: int) -> Chem.Mol:
        """
        Generates a polymer using a flexible list of monomer types and bead assignments.

        :param monomer_smiles: List of monomer SMILES.
        :param num_units: Corresponding number of monomer units per type.
        :return: Final polymer molecule.
        """
        monomer_iterators = []
        for monomer in self.monomer_smiles_list:
            monomer_result = self._create_monomer_residue(
                monomer
            )  # Get tuple (residue_mol, open_sites, bead_type)
            monomer_iterators.append(monomer_result)  # Create a cycling iterator

        monomer_iterator = cycle(monomer_iterators)

        first_residue, first_open_sites = next(monomer_iterator)

        polymer, idx = self._create_cap_residues(
            first_residue,
            open_sites=first_open_sites,
            use_open_site=0,
            add_to_map=True,
            add_to_sequence=True,
        )
        # print(self.debug_print_mol(polymer))
        for i in range(num_units - 2):
            residue_mol, open_sites = next(monomer_iterator)

            polymer, idx = self._add_monomer_to_polymer(
                polymer,
                monomer=residue_mol,
                prev_end_idx=idx,
                open_sites=open_sites,
                add_to_sequence=True,
            )

        # Add end monomer
        residue_mol, open_sites = next(monomer_iterator)
        end_monomer, _ = self._create_cap_residues(
            monomer=residue_mol,
            open_sites=open_sites,
            use_open_site=1,
            add_to_map=False,
            add_to_sequence=False,
        )
        polymer, _ = self._add_monomer_to_polymer(
            polymer,
            monomer=end_monomer,
            prev_end_idx=idx,
            open_sites=open_sites,
            add_to_sequence=True,
        )

        Chem.SanitizeMol(polymer)
        return polymer

    def _generate_filename(self, num_units: int) -> str:
        sanitised_smiles = [self.sanitize_filename(s) for s in self.monomer_smiles_list]
        file_name = "_".join(sanitised_smiles).lower()
        return f"{file_name}_{num_units}"

    @staticmethod
    def sanitize_filename(smiles: str) -> str:
        """
        Converts a SMILES string into a valid filename by replacing
        special characters with underscores.

        :param smiles: SMILES string to be sanitized
        :return: Sanitized filename string
        """
        return re.sub(r"[^\w\-_]", "_", smiles)  # Replace non-alphanumeric chars
