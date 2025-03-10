from modules.file_conversion.converter_factory import ConverterFactory
from modules.utils.shared.file_utils import (
    overwrite_directory,
    check_file_exists,
    check_directory_exists,
)
from modules.command_line_operation import CommandLineOperation
from config.acpype_config import (
    AcpypeOutputConfig,
)
from config.data_models.output_types import GromacsPaths
from modules.acpype.acpype_utils import (
    generate_acpype_paths,
    copy_acpype_files,
    rename_acpype_paths,
)
from config.paths import ACPYPE_POLYMER_NAME, TEMP_DIR
import os
from typing import Optional


class ACPYPEParameterizer(CommandLineOperation):

    def __init__(
        self,
        acpype_molecule_name: str = ACPYPE_POLYMER_NAME,
        temp_dir: str = TEMP_DIR,
    ):
        self.molecule_name = acpype_molecule_name
        self._temp_dir = temp_dir

        super().__init__()

    @property
    def temp_dir(self):
        return self._temp_dir

    @temp_dir.setter
    def temp_dir(self, value):
        self._temp_dir = value
        # Automatically update dependent attributes
        self.raw_output_dir = os.path.join(value, f"{self.molecule_name}.acpype")

    @property
    def raw_output_dir(self):
        return os.path.join(self.temp_dir, f"{self.molecule_name}.acpype")

    def acpype_command(self, input_file_path: str):
        acpype_command = [
            "acpype",
            "-i",
            input_file_path,
            "-o",
            "gmx",
            "-n",
            "0",
            "-a",
            "gaff2",
            "-b",
            self.molecule_name,
        ]

        return acpype_command

    def run(
        self,
        input_file_path: str,
        output_dir: str,
        acpype_output_config: AcpypeOutputConfig,
        new_file_name: Optional[str] = None,
        verbose: bool = False,
        skip_existing: bool = False,
    ) -> GromacsPaths:
        check_file_exists(input_file_path)
        check_directory_exists(output_dir)
        input_file_path = os.path.abspath(input_file_path)
        command = self.acpype_command(input_file_path=input_file_path)
        overwrite_directory(self.raw_output_dir)
        self._execute(command, cwd=self.temp_dir, verbose=verbose)

        raw_acpype_paths = generate_acpype_paths(
            acpype_output_config, self.raw_output_dir, self.molecule_name
        )
        copied_acpype_paths = copy_acpype_files(
            raw_acpype_paths, output_dir, skip_if_exists=skip_existing
        )
        if new_file_name:
            renamed_acpype_paths = rename_acpype_paths(
                copied_acpype_paths, new_file_name, suppress_warning=skip_existing
            )

            final_acpype_paths = renamed_acpype_paths
        else:
            final_acpype_paths = copied_acpype_paths

        return final_acpype_paths
