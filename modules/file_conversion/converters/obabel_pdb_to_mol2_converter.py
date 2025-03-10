from modules.file_conversion.converters.base_converter import (
    BaseConverter,
)
from modules.utils.shared.file_utils import (
    construct_output_file_path,
)
from typing import Optional, List, Tuple


class OBabelPDBtoMOL2Converter(BaseConverter):
    input_file_type = "pdb"
    output_file_type = "mol2"
    program = "Open Babel"

    def __init__(self):
        super().__init__()

    def _create_obabel_command(
        self, input_file_path: str, output_dir: Optional[str] = None
    ) -> Tuple[List, str]:
        output_file_path = construct_output_file_path(
            file_path=input_file_path,
            new_output_dir=output_dir,
            new_file_extension=self.output_file_type,
            suppress_warning=True,
            expected_file_type=self.input_file_type,
        )
        command = [
            "obabel",
            input_file_path,
            "-O",
            output_file_path,
            "--gen3D",
            "--h",  # add 3D coords via --h
        ]
        return command, output_file_path

    def _run_impl(
        self,
        input_file_path: str,
        output_dir: Optional[str] = None,
        verbose: bool = False,
    ):
        command, output_file_path = self._create_obabel_command(
            input_file_path, output_dir
        )
        self._execute(command, verbose=verbose)
        return output_file_path
