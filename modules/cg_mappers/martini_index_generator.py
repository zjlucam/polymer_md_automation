from modules.cg_mappers.base_map_generator import BaseMapGenerator
from typing import Optional


class MARTINIIndexGenerator(BaseMapGenerator):
    """
    Generates a MARTINI-compatible mapping file.
    """

    map_extension = "ndx"

    def __init__(self, polymer):
        super().__init__(polymer)

    def _generate_mapping(self, start_index: Optional[int] = None) -> str:
        """
        Generates the .ndx file content.
        """
        output = []
        bead_groups = {}

        # Collect all atom indices grouped by bead type
        for bead in self.bead_mappings:
            bead_type = bead["unique_name"]
            atom_indices = bead["atom_indices"]

            if bead_type not in bead_groups:
                bead_groups[bead_type] = []
            bead_groups[bead_type].extend(atom_indices)

        # Format each bead group into an .ndx section
        for bead_type, indices in bead_groups.items():
            output.append(
                f"[ {bead_type} ]\n" + " ".join(map(str, sorted(indices))) + "\n"
            )

        return "\n".join(output)
