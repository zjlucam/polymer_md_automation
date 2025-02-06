from typing import Dict, Tuple
from config.data_models.solvent import Solvent
from config.data_models.output_types import GromacsPaths
from modules.directory_parser.base_directory_traverser import BaseDirectoryParser
import os


class SolventDirectoryParser(BaseDirectoryParser):
    def parse_subdirectory(self, subdir_path: str) -> Tuple[GromacsPaths, Solvent]:

        itp_file = os.path.join(subdir_path, "SVT_GMX.itp")
        gro_file = os.path.join(subdir_path, "SVT_GMX.gro")
        top_file = os.path.join(subdir_path, "SVT_GMX.top")
        sol_properties_file = os.path.join(subdir_path, "sol_properties.txt")

        properties = self._parse_properties_file(sol_properties_file)

        required_fields = {
            "name",
            "compressibility",
            "molecular_weight",
            "density",
            "pdb_path",
            "SMILES",
        }
        if not required_fields.issubset(properties.keys()):
            print(
                f"Warning: Missing required fields in {sol_properties_file}. Skipping..."
            )
            return None

        try:
            compressibility = float(properties["compressibility"])
            density = float(properties["density"])
        except ValueError:
            print(
                f"Error: Invalid numeric values in {sol_properties_file}. Skipping..."
            )
            return None

        gromacs_paths = GromacsPaths(
            itp_path=itp_file if os.path.exists(itp_file) else None,
            gro_path=gro_file if os.path.exists(gro_file) else None,
            top_path=top_file if os.path.exists(top_file) else None,
        )

        solvent = Solvent(
            name=properties["name"],
            compressibility=compressibility,
            density=density,
            pdb_path=properties["pdb_path"],
            molecular_weight=(float(properties["molecular_weight"])),
        )

        solvent_smiles = properties["SMILES"]

        return gromacs_paths, solvent, solvent_smiles

    @staticmethod
    def _parse_properties_file(file_path: str) -> Dict[str, str]:

        properties = {}
        if os.path.exists(file_path):
            with open(file_path, "r") as file:
                for line in file:
                    line = line.strip()
                    if "=" in line:
                        key, value = line.split("=", 1)
                        properties[key.strip()] = value.strip()
        return properties
