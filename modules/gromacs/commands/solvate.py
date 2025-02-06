from modules.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from modules.utils.shared.file_utils import (
    check_file_type,
    check_directory_exists,
)
from typing import Optional, List, Dict
import os


class Solvate(BaseGromacsCommand):

    def __init__(self):
        super().__init__()

    def run(
        self,
        solute_gro_path: str,
        solvent_gro_path: str,
        input_topol_path: str,
        output_dir: str,
        output_name: str = "solvated.gro",
    ):
        check_file_type(solute_gro_path, "gro")
        check_file_type(solvent_gro_path, "gro")
        check_file_type(input_topol_path, "top")
        check_directory_exists(output_dir)
        output_gro_path = os.path.join(output_dir, output_name)

        command = self._create_command(
            solute_gro_path=solute_gro_path,
            solvent_gro_path=solvent_gro_path,
            input_topol_path=input_topol_path,
            output_gro_path=output_gro_path,
        )

        self._execute(command)

        return output_gro_path

    def _create_command(
        self,
        solute_gro_path: str,
        solvent_gro_path: str,
        input_topol_path: str,
        output_gro_path: str,
    ) -> List[str]:

        solvate_command = [
            "gmx",
            "solvate",
            "-cp",
            solute_gro_path,
            "-cs",
            solvent_gro_path,
            "-p",
            input_topol_path,
            "-o",
            output_gro_path,
        ]
        return solvate_command
