from modules.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from modules.utils.shared.file_utils import (
    check_file_type,
    prepare_output_file_path,
)
from config.constants import MassUnits
from typing import Optional, Dict, Tuple, List


class InsertMolecules(BaseGromacsCommand):
    default_mass_units = MassUnits.GRAM

    def __init__(self):
        super().__init__()

    def run(
        self,
        input_box_gro_path: str,
        solvent_gro_path: str,
        num_molecules: int,
        output_dir: Optional[str] = None,
        output_name: Optional[str] = None,
    ):
        check_file_type(input_box_gro_path, "gro")
        check_file_type(solvent_gro_path, "gro")
        output_path = prepare_output_file_path(
            input_box_gro_path,
            output_extension="gro",
            output_dir=output_dir,
            output_name=output_name,
        )
        command, output_path = self._create_command(
            input_box_gro_path=input_box_gro_path,
            solvent_gro_path=solvent_gro_path,
            num_molecules=num_molecules,
            output_path=output_path,
        )

        self._execute(command)
        return output_path

    def _create_command(
        self,
        input_box_gro_path: str,
        solvent_gro_path: str,
        num_molecules: int,
        output_path: str = None,
    ) -> Tuple[List[str], str]:

        if not output_path:
            output_path = input_box_gro_path
            append_flag = "-append"
        else:
            append_flag = None

        # Build the GROMACS command
        command = [
            "gmx",
            "insert-molecules",
            "-f",
            input_box_gro_path,  # Input structure file with the existing solvent
            "-ci",
            solvent_gro_path,  # Input solvent file to insert
            "-nmol",
            str(num_molecules),  # Number of molecules to insert
            "-o",
            output_path,  # Output file
        ]

        # Add the append flag if applicable
        if append_flag:
            command.append(append_flag)

        return command, output_path
