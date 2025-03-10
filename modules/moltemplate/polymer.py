from modules.moltemplate.base_molecule import BaseMoltemplateMolecule
from typing import Dict, Optional
from modules.lammps.parsers.open_mscg_data_parser import OpenMSCGDataParser
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MoltemplatePolymer(BaseMoltemplateMolecule):

    def __init__(
        self,
        n_units: int,
        polymer_generator: BasePolymerGenerator,
        openmscg_parser: OpenMSCGDataParser,
        replace_using_openmscg_mass_map: bool = False,
    ):
        super().__init__(
            openmscg_parser=openmscg_parser,
            use_openmscg_mass_map=replace_using_openmscg_mass_map,
        )
        self.n_units = n_units
        self.polymer_generator = polymer_generator
        self.terminal_left = None
        self.terminal_right = None
        self.middle_beads = None
        self.cg_map = polymer_generator.cg_map
        self.n_repeats = None
        self.actual_n = None
        self.openmscg_parser = openmscg_parser
        self.unique_bead_types = None
        self.lt_block = self.generate_lt()

    def _get_mass_mapping(self) -> None:
        unique_bead_types = []
        mass_mapping = {}
        for bead in self.cg_map:
            bead_type = bead["bead_type"]
            mass = sum(bead["x-weight"])

            if bead_type not in unique_bead_types:
                unique_bead_types.append(bead_type)

            if bead_type not in self.mass_mapping:
                mass_mapping[bead_type] = mass
            else:
                stored_mass = mass_mapping[bead_type]

                if stored_mass == 0 and mass != 0:
                    mass_mapping[bead_type] = mass

                elif stored_mass != 0 and mass != 0 and abs(mass - stored_mass) > 1e-12:
                    logger.warning(
                        f"Bead type '{bead_type}' has conflicting masses: "
                        f"{stored_mass} vs {mass}. Retaining {stored_mass}"
                    )

        self.terminal_left = unique_bead_types[0]
        self.terminal_right = unique_bead_types[-1]
        self.middle_beads = unique_bead_types[1:-1]

        self.unique_bead_types = unique_bead_types
        self.mass_mapping = mass_mapping
        return mass_mapping

    def _round_polymer_length(self) -> int:
        repeating_units = len(self.middle_beads)
        n_repeats = (self.n_units - 2) // repeating_units
        actual_n = n_repeats * repeating_units + 2
        logging.info(f"Desired number of beads: {self.n_units}")
        logging.info("Rounding...")
        logging.info(f"Creating polymer with beads: {actual_n}")
        return n_repeats, actual_n

    def _retrieve_bond_length(self, bead_1: str, bead_2: str) -> float:
        return self.openmscg_parser.retrieve_bond_length(bead_1=bead_1, bead_2=bead_2)

    def _generate_lt(self) -> str:
        lines = ["Molecule Polymer {"]
        x_pos = 0.0

        lines.append(f"  $atom:{self.terminal_left} {x_pos:.1f} 0.0 0.0")
        last_bead = self.terminal_left

        self.n_repeats, self.actual_n = self._round_polymer_length()

        for i in range(self.n_repeats):
            for bead in self.middle_beads:
                x_pos += self._retrieve_bond_length(bead_1=last_bead, bead_2=bead)
                lines.append(f"  $atom:{bead}_{i+1} {x_pos:.1f} 0.0 0.0")
                last_bead = bead

        x_pos += self._retrieve_bond_length(
            bead_1=last_bead, bead_2=self.terminal_right
        )
        lines.append(f"  $atom:{self.terminal_right} {x_pos:.1f} 0.0 0.0")

        lines.append("\n  # Bonds")
        all_atoms = (
            [self.terminal_left]
            + [f"{b}_{i+1}" for i in range(self.n_repeats) for b in self.middle_beads]
            + [self.terminal_right]
        )
        for i in range(len(all_atoms) - 1):
            lines.append(
                f"  bond:bond_{i+1} @atom:{all_atoms[i]} @atom:{all_atoms[i+1]}"
            )

        lines.append("\n  # Angles")
        for i in range(len(all_atoms) - 2):
            lines.append(
                f"  angle:angle_{i+1} @atom:{all_atoms[i]} @atom:{all_atoms[i+1]} @atom:{all_atoms[i+2]}"
            )

        lines.append("}\n")
        return "\n".join(lines)
