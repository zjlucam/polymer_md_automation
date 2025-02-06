import os

# from A_config.paths import GROMACS_OUTPUT_SUBDIR
from config.constants import DEFAULT_BOX_SIZE_NM
from typing import Optional, List
from modules.utils.shared.calculation_utils import calculate_num_particles
from modules.utils.shared.file_utils import (
    check_directory_exists,
    check_file_does_not_exist,
)
from modules.gromacs.commands.base_gromacs_command import BaseGromacsCommand
import logging
from typing import Tuple
from config.constants import MassUnits
import subprocess
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.gro_handler import GroHandler

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class SolventInsertion(BaseGromacsCommand):
    output_name = "solvent_box.gro"
    default_mass_units = MassUnits.GRAM

    def __init__(self):
        super().__init__()

    def run(
        self,
        input_pdb_path: str,
        output_dir: str,
        desired_density: float,
        molecular_weight: float,
        box_size_nm: List[float] = DEFAULT_BOX_SIZE_NM,
        num_iterations_max: int = 10,
        tolerance: float = 0.05,
        verbose: bool = False,
    ) -> str:
        check_directory_exists(output_dir)
        self._validate_parameters(
            desired_density, molecular_weight, num_iterations_max, tolerance
        )

        gro_path = os.path.join(output_dir, self.output_name)
        check_file_does_not_exist(gro_path, suppress_error=True, delete_file=True)

        target_num_particles = self._calculate_target_particles(
            box_size_nm, molecular_weight, desired_density, verbose
        )

        for iteration in range(1, num_iterations_max + 1):
            if self._is_target_achieved(
                gro_path, target_num_particles, tolerance, iteration, verbose
            ):

                return gro_path

            remaining_particles = self._calculate_remaining_particles(
                gro_path, target_num_particles
            )
            self._add_particles(
                input_pdb_path, gro_path, box_size_nm, remaining_particles, verbose
            )

        # Raise an error if target is not met after max iterations
        raise RuntimeError(
            f"Failed to reach the target number of particles after {num_iterations_max} iterations."
        )

    def _calculate_target_particles(
        self, box_size_nm, molecular_weight, desired_density, verbose
    ):
        """
        Calculates the target number of particles based on box size, molecular weight, and desired density.
        """
        target_num_particles = calculate_num_particles(
            box_size_nm,
            molecular_weight=molecular_weight,
            density_SI=desired_density,
            box_units=self.default_units,
            mass_units=self.default_mass_units,
        )
        if verbose:
            logger.info(f"Target number of particles: {round(target_num_particles)}")
        return target_num_particles

    def _is_target_achieved(
        self, gro_path, target_num_particles, tolerance, iteration, verbose
    ):
        """
        Checks if the target number of particles has been achieved.
        """
        current_num_particles = (
            self._count_particles(gro_path) if os.path.exists(gro_path) else 0
        )
        remaining_particles = target_num_particles - current_num_particles

        if verbose:
            logger.info(
                f"Iteration {iteration}: Current particles: {current_num_particles}, Remaining: {remaining_particles}"
            )

        return remaining_particles <= target_num_particles * tolerance

    def _calculate_remaining_particles(self, gro_path, target_num_particles):
        """
        Calculates the number of remaining particles needed to reach the target.
        """
        current_num_particles = (
            self._count_particles(gro_path) if os.path.exists(gro_path) else 0
        )
        return target_num_particles - current_num_particles

    def _add_particles(
        self, input_pdb_path, gro_path, box_size_nm, remaining_particles, verbose
    ):
        """
        Adds remaining particles to the system by executing the insert command.
        """
        # Check if the output file exists
        append = os.path.exists(gro_path)

        # Create the command with or without the append flag
        command, _ = self._create_insert_command(
            input_pdb_path=input_pdb_path,
            box_size_nm=box_size_nm,
            num_particles=remaining_particles,
            output_path=gro_path,
            append=append,  # Only append if the file exists
        )
        self._execute(command, verbose=verbose)

    def _create_insert_command(
        self,
        input_pdb_path: str,
        box_size_nm: List[float],
        num_particles: int,
        output_path: str,
        append: bool = False,
    ) -> Tuple[List[str], str]:
        command = [
            "gmx",
            "insert-molecules",
            "-ci",
            input_pdb_path,
            "-nmol",
            str(num_particles),
            "-box",
            str(box_size_nm[0]),
            str(box_size_nm[1]),
            str(box_size_nm[2]),
            "-o",
            output_path,
        ]
        if append:
            command.append("-f")
            command.append(output_path)  # Append to the existing file
        return command, output_path

    def _validate_parameters(
        self,
        desired_density: float,
        molecular_weight: float,
        num_iterations_max: int,
        tolerance: float,
    ):
        if desired_density <= 0:
            raise ValueError("Desired density must be positive.")
        if molecular_weight <= 0:
            raise ValueError("Molecular weight must be positive.")
        if num_iterations_max <= 0:
            raise ValueError("Number of iterations must be positive.")
        if not (0 < tolerance < 1):
            raise ValueError("Tolerance must be between 0 and 1.")

    def _count_particles(self, gro_file: str) -> int:
        """
        Count the number of molecules (particles) in a .gro file.

        Args:
            gro_file (str): Path to the `.gro` file.

        Returns:
            int: Number of unique molecules (residues) in the file.
        """
        if not os.path.exists(gro_file):
            return 0

        # Use FileSplitter and GroHandler to parse the .gro file
        parser = GromacsParser()
        sections = parser.parse(gro_file)

        # Assuming the `.gro` file content is the first section
        gro_section = next(iter(sections.values()))
        gro_handler = GroHandler()
        gro_handler.process(gro_section)

        residue_numbers = gro_handler.content["Residue Number"].unique()
        return len(residue_numbers)
