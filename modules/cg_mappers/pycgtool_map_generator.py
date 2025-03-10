from modules.cg_mappers.base_map_generator import BaseMapGenerator
import re
from typing import Optional
from modules.gromacs.parsers.itp_parser import ITPParser


class PyCGToolMapGenerator(BaseMapGenerator):
    """
    Generates a MARTINI-compatible mapping file.
    """

    map_extension = "map"

    def __init__(
        self,
        polymer,
        molecule_name: Optional[str] = "UNL",
    ):
        super().__init__(polymer)
        if not molecule_name:
            self._molecule_name = polymer.res_name
        else:
            self._molecule_name = molecule_name
        self._polymer_content = None
        self._solvent_content = None

    @property
    def molecule_name(self) -> str:
        return self._molecule_name

    @molecule_name.setter
    def molecule_name(self, name: str):
        self._molecule_name = name

    def _generate_polymer_mapping(self, start_index: Optional[int] = None) -> str:
        """
        Generates the MARTINI mapping file content.
        """
        output = [f"[{self._molecule_name}]"]

        atom_lines = []
        atom_index = 1

        for bead in self.bead_mappings:
            bead_name = bead["unique_name"]
            bead_type = bead["bead_type"]

            atom_names = [
                self.reformat_atom_name(atom, start_index=start_index)
                for atom in bead["atom_names"]
            ]

            # **Group all atom names per bead in a single line**
            atom_line = f"{bead_name}\t{bead_type}\t" + " ".join(atom_names)
            output.append(atom_line)

        self._polymer_content = "\n".join(output) + "\n"

    def add_solvent_to_map(
        self, itp_file: str, bead_name: str = "SOL", start_index: Optional[int] = None
    ) -> str:
        """
        Generates a coarse-grained mapping from a solvent .itp file.
        Each solvent molecule is mapped to one bead.

        Args:
            itp_file (str): Path to the solvent .itp file.
            bead_name (str): The CG bead name to assign.

        Returns:
            str: The formatted CG mapping.
        """
        atom_section = False
        solvent_name = None
        mapping_lines = []

        itp_data = ITPParser(itp_file)
        itp_df = itp_data.retrieve_content("atoms")
        solvent_name = itp_df["res"].iloc[0]  # Take the first occurrence

        atom_names = " ".join(
            itp_df["atom"]
            .apply(lambda x: self.reformat_atom_name(x, start_index))
            .tolist()
        )

        output = f"[ {solvent_name} ]\n{bead_name}\t{bead_name}\t{atom_names}\n"

        self._solvent_content = output

    def _generate_mapping(self, start_index: Optional[str] = None) -> str:
        """
        Generates the MARTINI mapping file content.
        """
        self._generate_polymer_mapping(start_index=start_index)

        if self._solvent_content:
            return self._polymer_content + self._solvent_content

        return self._polymer_content
