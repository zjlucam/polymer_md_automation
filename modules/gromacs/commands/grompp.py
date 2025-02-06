from modules.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from modules.utils.shared.file_utils import (
    check_file_type,
    prepare_output_file_path,
)
from typing import Optional, List


class Grompp(BaseGromacsCommand):

    def __init__(self):
        super().__init__()

    def run(
        self,
        mdp_file_path: str,
        input_gro_path: str,
        input_topol_path: str,
        output_dir: Optional[str] = None,
        output_name: Optional[str] = None,
        verbose: bool = False,
    ) -> str:
        output_tpr_path = prepare_output_file_path(
            input_gro_path,
            output_extension="tpr",
            output_dir=output_dir,
            output_name=output_name,
        )
        command = self._create_command(
            mdp_file_path,
            input_gro_path,
            input_topol_path,
            output_tpr_path,
        )
        self._execute(command, verbose=verbose)
        return output_tpr_path

    def _create_command(
        self,
        mdp_file_path: str,
        input_gro_path: str,
        input_topol_path: str,
        output_tpr_path: str,
    ) -> List[str]:
        check_file_type(mdp_file_path, "mdp")
        check_file_type(input_gro_path, "gro")
        check_file_type(input_topol_path, "top")
        check_file_type(output_tpr_path, "tpr")

        return [
            "gmx",
            "grompp",
            "-f",
            mdp_file_path,
            "-c",
            input_gro_path,
            "-p",
            input_topol_path,
            "-o",
            output_tpr_path,
        ]
