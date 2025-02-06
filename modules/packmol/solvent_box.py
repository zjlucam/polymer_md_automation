from typing import Dict, List
from modules.packmol.base_packmol_operation import BasePackmolOperation
from config.data_models.solvent import Solvent
from modules.utils.shared.calculation_utils import (
    calculate_num_particles,
    box_dimensions_check_wrapper,
)
from config.constants import LengthUnits, MassUnits
from modules.utils.shared.file_utils import file_type_check_wrapper
<<<<<<< HEAD
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
=======
>>>>>>> 91758eb (cleaned up)


class PackmolSolventBox(BasePackmolOperation):
    template_name = "solvent_box_template.inp"

    def generate_template_params(
        self,
        input_pdb: str,
        output_file: str,
        solvent: Solvent,
        box_size_nm: List[float],
        **kwargs,
    ) -> Dict[str, str]:

        box_size_angstrom = [
            dim * LengthUnits.NANOMETER.value / LengthUnits.ANGSTROM.value
            for dim in box_size_nm
        ]

<<<<<<< HEAD
        logger.info(f"Box_size_angstrom: {box_size_angstrom}")
=======
        print("box_size_angstrom", box_size_angstrom)
>>>>>>> 91758eb (cleaned up)
        num_molecules = calculate_num_particles(
            box_dimensions=box_size_angstrom,
            molecular_weight=solvent.molecular_weight,
            density_SI=solvent.density,
            box_units=LengthUnits.ANGSTROM,
            mass_units=MassUnits.GRAM,
        )
<<<<<<< HEAD
        logger.info(f"Number of molecules added: {num_molecules}")
=======
        print(num_molecules)
>>>>>>> 91758eb (cleaned up)
        return {
            "output_file": output_file,
            "solvent_file": input_pdb,
            "num_molecules": num_molecules,
            "box_size_x": box_size_angstrom[0],
            "box_size_y": box_size_angstrom[1],
            "box_size_z": box_size_angstrom[2],
        }
