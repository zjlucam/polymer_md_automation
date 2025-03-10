from modules.command_line_operation import CommandLineOperation

from abc import ABC, abstractmethod
from typing import List, Dict, Optional
import os
from config.paths import TEMP_DIR
import logging
from modules.utils.shared.file_utils import (
    file_type_check_wrapper,
    check_file_exists,
    directory_exists_check_wrapper,
)
from config.paths import PACKMOL_TEMPLATE_DIR
import subprocess

logger = logging.getLogger(__name__)


class BasePackmolOperation(CommandLineOperation, ABC):
    """
    Abstract base class for operations using Packmol. Provides a framework
    for generating Packmol input scripts and running Packmol commands.
    """

    template_name: str
    TEMPLATE_DIR = PACKMOL_TEMPLATE_DIR

    def __init_subclass__(cls):
        super().__init_subclass__()
        if not hasattr(cls, "template_name"):
            raise TypeError(f"Class {cls.__name__} must define 'template_name'.")

    def __init__(
        self,
        tolerance: float = 2.0,
        filetype: str = "pdb",
    ):
        super().__init__()
        self.tolerance = tolerance
        self.filetype = filetype

    @abstractmethod
    def generate_template_params(
        self, input_pdb: str, output_file: str, **kwargs
    ) -> Dict[str, str]:
        """
        Abstract method for subclasses to define template parameters.
        """
        pass

    def generate_input_script(
        self,
        template_params: Dict[str, str],
        temp_output_dir: str = TEMP_DIR,
        output_name: str = "packmol_input",
        template_path: Optional[str] = None,
    ) -> str:
        """
        Generate the Packmol input script from a template.

        Args:
            template_params (Dict[str, str]): Parameters to populate the template.
            temp_output_dir (str): Directory to store the temporary input script.
            output_name (str): Name of the output file (without extension).
            template_path (Optional[str]): Path to a custom template file.

        Returns:
            str: Path to the generated input script.
        """
        # Ensure `temp_output_dir` exists
        os.makedirs(temp_output_dir, exist_ok=True)

        # Determine the template file
        if not template_path:
            template_path = os.path.join(self.TEMPLATE_DIR, self.template_name)

        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Template file not found: {template_path}")

        # Read and populate the template
        with open(template_path, "r") as f:
            template_content = f.read()

        try:
            input_script = template_content.format(**template_params)
        except KeyError as e:
            raise ValueError(f"Missing required template parameter: {e}")

        # Save the script to the specified temp output directory
        script_path = os.path.join(temp_output_dir, f"{output_name}.inp")
        with open(script_path, "w", newline="\n") as f:
            f.write(input_script)

        logger.info(f"Generated Packmol input script at {script_path}")
        return script_path

    def run(
        self,
        input_pdb: str,
        output_name: Optional[str] = None,
        output_dir: str = TEMP_DIR,
        temp_output_dir: str = TEMP_DIR,
        template_path: Optional[str] = None,
        verbose: bool = False,
        **kwargs,
    ) -> str:
        if not output_name:
            output_name = (
                f"{os.path.splitext(os.path.basename(input_pdb))[0]}_packed.pdb"
            )
        logger.info(f"Starting Packmol operation: {self.__class__.__name__}")

        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Step 1: Generate template parameters
        template_params = self.generate_template_params(
            input_pdb=input_pdb,
            output_file=os.path.join(output_dir, output_name),
            **kwargs,
        )
        logger.debug(f"Generated template parameters: {template_params}")

        # Step 2: Generate the Packmol input script
        input_script = self.generate_input_script(
            template_params=template_params,
            temp_output_dir=temp_output_dir,
            output_name=os.path.splitext(os.path.basename(output_name))[0],
            template_path=template_path,
        )
        logger.debug(f"Generated Packmol input script: {input_script}")

        # Step 3: Run Packmol with the generated input script
        command = self._create_packmol_command(input_script)

        ##################################3
        # self._execute(command, verbose=verbose)
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

        #############################################
        # Step 4: Validate output file existence
        final_output_path = os.path.join(output_dir, output_name)
        if not os.path.exists(final_output_path):
            logger.error(
                f"Packmol failed to generate the output file: {final_output_path}"
            )
            raise RuntimeError(
                f"Output file not found after Packmol execution: {final_output_path}"
            )

        logger.info(
            f"Packmol operation completed successfully. Output file: {final_output_path}"
        )
        return final_output_path

    def _create_packmol_command(self, input_script: str) -> List[str]:
        """
        Create the command to run Packmol.

        Args:
            input_script (str): Path to the input script.

        Returns:
            List[str]: Command to execute Packmol.
        """
        return f"packmol < {input_script}"
