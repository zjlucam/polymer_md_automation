import os
import subprocess
from modules.file_conversion.converters.base_converter import BaseConverter
from typing import Optional, List, Tuple
from modules.gromacs.commands.editconf import Editconf
from pathlib import Path


class EditconfGROtoPDBConverter(BaseConverter):
    input_file_type = "gro"
    output_file_type = "pdb"
    program = "GROMACS"

    def __init__(self):
        super().__init__()

    def _run_impl(
        self,
        input_file_path: str,
        output_dir: Optional[str] = None,
        verbose: bool = False,
        box_size_nm: Optional[List[float]] = None,
        output_name: Optional[str] = None,
    ):
        editconf = Editconf()
        if output_name is None:
            output_name = Path(input_file_path).stem + ".pdb"
        output_path = editconf.run(
            input_file_path, output_dir, box_size_nm, output_name
        )
        return output_path
