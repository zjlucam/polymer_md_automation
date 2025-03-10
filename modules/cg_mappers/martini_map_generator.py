from modules.cg_mappers.base_map_generator import BaseMapGenerator
from typing import Optional


class MARTINIMapGenerator(BaseMapGenerator):
    """
    Generates a MARTINI-compatible mapping file.
    """

    map_extension = "map"

    def __init__(self, polymer):
        super().__init__(polymer)

    def _generate_mapping(self, start_index: Optional[str] = None) -> str:
        """
        Generates the MARTINI mapping file content.
        """
        output = ["[ to ]", "martini\n"]

        # Extract unique bead names from cg_map
        bead_types = sorted(set(bead["unique_name"] for bead in self.bead_mappings))
        output.append("[ martini ]\n" + " ".join(bead_types) + "\n\n")

        # Generate [ atoms ] section
        output.append("[ atoms ]")
        atom_lines = []
        atom_index = 1

        for bead in self.bead_mappings:
            bead_type = bead["unique_name"]
            atom_indices = bead["atom_indices"]
            atom_names = bead["atom_names"]

            for idx, atom_name in zip(atom_indices, atom_names):
                atom_name = self.reformat_atom_name(atom_name)
                atom_lines.append(f"{atom_index}\t{atom_name}\t{bead_type}")
                atom_index += 1

        output.append("\n".join(atom_lines))
        return "\n".join(output)
