from modules.rdkit.base_molecule_generator import BaseMoleculeGenerator
from config.paths import SOLVENT_PDB_DIR, TEMP_DIR
from modules.acpype.acpype_parametizer import ACPYPEParameterizer
from modules.utils.shared.file_utils import check_directory_exists
from modules.file_conversion.converters.obabel_pdb_to_mol2_converter import (
    OBabelPDBtoMOL2Converter,
)
from typing import Optional, List
import os
import logging
from modules.rdkit.solvent_generator import SolventGenerator
from config.acpype_config import AcpypeOutputConfig

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SolventParametiser:
    sol_text_file = "sol_properties.txt"

    def __init__(
        self,
        solvent_name: str,
        solvent_smiles: str,
        solvent_compressibility: str,
        solvent_density: float,
        output_dir: str,
        verbose: bool = False,
        solvent_resname: str = "SOL",
    ):
        self.output_dir = os.path.join(output_dir, solvent_name)
        self.parametizer = ACPYPEParameterizer(acpype_molecule_name="SVT")
        self.verbose = verbose
        self.file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
        self.generator = SolventGenerator(
            solvent_name=solvent_name,
            solvent_smiles=solvent_smiles,
            solvent_compressibility=solvent_compressibility,
            solvent_density=solvent_density,
            solvent_resname=solvent_resname,
        )
        self.solvent = self.generator.generate_solvent(output_dir=self.output_dir)
        self.params = {
            "name": solvent_name,
            "compressibility": solvent_compressibility,
            "SMILES": solvent_smiles,
            "density": solvent_density,
            "pdb_path": self.solvent.pdb_path,
            "molecular_weight": self.solvent.molecular_weight,
        }

    def run(self, mol2_output_dir: str = TEMP_DIR):
        check_directory_exists(self.output_dir)
        check_directory_exists(mol2_output_dir)
        converter = OBabelPDBtoMOL2Converter()
        mol2_file = converter.run(
            self.solvent.pdb_path, mol2_output_dir, verbose=self.verbose
        )
        self.parametizer.run(
            mol2_file, self.output_dir, self.file_config, verbose=self.verbose
        )
        self.write_params_and_list_to_file(
            output_dir=self.output_dir,
            text_file_name=self.sol_text_file,
            params=self.params,
        )

    def write_params_and_list_to_file(
        self,
        output_dir: str,
        text_file_name: str,
        params: dict,
        content_list: Optional[List[str]] = None,
    ):

        file_path = os.path.join(output_dir, text_file_name)

        os.makedirs(output_dir, exist_ok=True)

        with open(file_path, "w") as file:

            for key, value in params.items():
                file.write(f"{key}={value}\n")

            if content_list:
                file.write("---LIST---\n")

                # Write the list items
                file.writelines(f"{item}\n" for item in content_list)
