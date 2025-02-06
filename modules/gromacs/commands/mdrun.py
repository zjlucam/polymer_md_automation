from modules.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from typing import Optional, Dict, Tuple, List
import os


# NOTE: NEED TO FIX, think abt temp, and where to input/output dir - e.g. in the new base class for equilibrium or not
class MDrun(BaseGromacsCommand):
    def __init__(self):
        super().__init__()

    def run(
        self,
        input_tpr_path: str,
        output_name: str,
        verbose: bool = False,
        additional_flags: Optional[List[str]] = None,
    ) -> Dict[str, str]:
        command = self._create_command(input_tpr_path, output_name, additional_flags)
        self._execute(command, verbose=verbose)
        output_files = {
            ext: f"{output_name}.{ext}" for ext in ["gro", "log", "edr", "trr"]
        }

        return {k: v for k, v in output_files.items() if os.path.exists(v)}

    def _create_command(
        self,
        input_tpr_path: str,
        output_name: str,
        additional_flags: Optional[List[str]] = None,
    ) -> List[str]:
        command = ["gmx", "mdrun", "-s", input_tpr_path, "-deffnm", output_name]
        if additional_flags:
            command.extend(additional_flags)
        return command
