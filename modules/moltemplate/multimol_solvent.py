from modules.moltemplate.base_molecule import BaseMoltemplateMolecule
from modules.lammps.parsers.open_mscg_data_parser import OpenMSCGDataParser
from config.data_models.solvent import Solvent
from typing import Dict


class MultimolMoltemplateSolvent(BaseMoltemplateMolecule):

    def __init__(
        self,
        solvent: Solvent,
        openmscg_parser: OpenMSCGDataParser,
        sol_resname="SOL",
        replace_using_openmscg_mass_map: bool = True,
        sol_to_bead_ratio: int = 4,
    ):
        super().__init__(
            openmscg_parser=openmscg_parser,
            use_openmscg_mass_map=replace_using_openmscg_mass_map,
        )
        self.solvent = solvent
        self.sol_to_bead_ratio = sol_to_bead_ratio
        self.resname = sol_resname
        self.mass = self.solvent.molecular_weight * sol_to_bead_ratio
        self.mass_mapping[self.resname] = self.mass
        self.lt_block = self.generate_lt()

    def _get_mass_mapping(self) -> Dict[str, float]:
        return {self.resname: self.mass}

    def _generate_lt(self) -> str:
        return f"Molecule Solvent {{\n  $atom:{self.resname} 0.0 0.0 0.0\n}}\n"
