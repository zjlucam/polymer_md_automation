from modules.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from modules.utils.shared.file_utils import (
    check_file_type,
    prepare_output_file_path,
)
from typing import Optional, Dict, Tuple, List
import subprocess
from typing import Optional, List, Tuple
from modules.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from modules.utils.shared.file_utils import (
    check_file_type,
    prepare_output_file_path,
)


class GenIon(BaseGromacsCommand):
    def __init__(self):
        super().__init__()

    def _create_command(
        self,
        tpr_path: str,
        top_path: str,
        output_path: str,
        pname: str,
        nname: str,
        sol_name: str,  # Pass solvent name explicitly
        conc: Optional[float],
    ) -> Tuple[List[str], str]:

        command = [
            "gmx",
            "genion",
            "-s",
            tpr_path,
            "-o",
            output_path,
            "-p",
            top_path,
            "-pname",
            pname,
            "-nname",
            nname,
            "-neutral",
        ]

        if conc:
            command += ["-conc", str(conc)]

        return command, output_path

    def _generate_index_file(self, gro_file: str) -> str:
        """
        Generate an index file to ensure solvent selection is automated.
        """
        index_file = gro_file.replace(".gro", ".ndx")
        command = ["gmx", "make_ndx", "-f", gro_file, "-o", index_file]

        try:
            subprocess.run(command, input="q\n", text=True, check=True)
            return index_file
        except subprocess.CalledProcessError:
            raise RuntimeError("Failed to generate index file for solvent selection.")

    def run(
        self,
        input_box_gro_path: str,
        tpr_path: str,
        top_path: str,
        sol_name: str = "SOL",  # Default to SOL
        output_dir: Optional[str] = None,
        output_name: str = "neutralized_polymer",
        pname: str = "NA",
        nname: str = "CL",
        conc: Optional[float] = None,
        verbose: bool = True,
    ) -> str:
        """
        Runs genion while ensuring correct solvent selection.
        """
        check_file_type(input_box_gro_path, "gro")
        check_file_type(tpr_path, "tpr")
        check_file_type(top_path, "top")

        output_path = prepare_output_file_path(
            input_box_gro_path,
            "gro",
            output_dir,
            output_name,
        )

        # 🔄 Generate an index file to avoid manual input
        index_file = self._generate_index_file(input_box_gro_path)

        # 🔄 Run genion with explicit SOL group
        command, output_path = self._create_command(
            tpr_path=tpr_path,
            top_path=top_path,
            output_path=output_path,
            pname=pname,
            nname=nname,
            sol_name=sol_name,
            conc=conc,
        )

        # Ensure we select the solvent group automatically
        command += ["-n", index_file]
        command_input = f"{sol_name}\n"

        self._execute(command, input=command_input, verbose=verbose)
        return output_path
