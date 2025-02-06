from modules.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from modules.utils.shared.file_utils import check_file_exists, check_directory_exists
from typing import Optional, List, Dict
import os


class Editconf(BaseGromacsCommand):
    def __init__(self):
        super().__init__()

    def run(
        self,
        input_gro_or_pdb_path: str,
        output_dir: str,
        box_size_nm: Optional[List[float]] = None,
        output_name: str = "edited_box.gro",
    ):
        check_file_exists(input_gro_or_pdb_path)
        check_directory_exists(output_dir)
        output_gro_path = os.path.join(output_dir, output_name)

        command = self._create_command(
            input_gro_or_pdb_path=input_gro_or_pdb_path,
            output_gro_path=output_gro_path,
            box_size_nm=box_size_nm,
        )

        self._execute(command)
        return output_gro_path

    def _create_command(
        self,
        input_gro_or_pdb_path: str,
        output_gro_path: str,
        box_size_nm: Optional[List[float]],
    ) -> List[str]:
        command = [
            "gmx",
            "editconf",
            "-f",
            input_gro_or_pdb_path,
            "-o",
            output_gro_path,
        ]
        if box_size_nm is not None:
            command.extend(
                [
                    "-box",
                    str(box_size_nm[0]),
                    str(box_size_nm[1]),
                    str(box_size_nm[2]),
                ]
            )
        return command
