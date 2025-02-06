from config.data_models.output_types import GromacsPaths
from modules.directory_parser.base_directory_traverser import BaseDirectoryParser
from typing import List, Tuple, Optional
import os


class PolymerDirectoryParser(BaseDirectoryParser):

    def parse_subdirectory(
        self, subdir_path: str
    ) -> Tuple[GromacsPaths, List[str], int]:
        # Define expected file paths
        itp_file = os.path.join(subdir_path, "POLY_GMX.itp")
        gro_file = os.path.join(subdir_path, "POLY_GMX.gro")
        top_file = os.path.join(subdir_path, "POLY_GMX.top")
        smiles_file = os.path.join(subdir_path, "smiles.txt")

        # Read SMILES and num_units
        monomer_smiles_list, num_units = self._parse_smiles_file(smiles_file)

        if num_units is None or not monomer_smiles_list:
            print(
                f"Warning: Skipping {subdir_path} due to missing or malformed smiles.txt."
            )
            return None

        # Create GromacsPaths object
        gromacs_paths = GromacsPaths(
            itp_path=itp_file if os.path.exists(itp_file) else None,
            gro_path=gro_file if os.path.exists(gro_file) else None,
            top_path=top_file if os.path.exists(top_file) else None,
        )

        return gromacs_paths, monomer_smiles_list, num_units

    @staticmethod
    def _parse_smiles_file(file_path: str) -> Tuple[List[str], Optional[int]]:

        if not os.path.exists(file_path):
            return [], None

        with open(file_path, "r") as file:
            lines = [line.strip() for line in file if line.strip()]
            if not lines:
                return [], None
            try:
                num_units = int(lines[0])
                monomer_smiles_list = lines[1:]
                return monomer_smiles_list, num_units
            except ValueError:
                return [], None
