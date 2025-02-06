from modules.rdkit.solvent_generator import SolventGenerator
from modules.workflows.atomistic.no_parameterising_equilibriator import (
    NoParametiserPolymerEquilibriationWorkflow,
)
from modules.gromacs.analyser import GromacsAnalyser
from typing import List, Any, Dict
from modules.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)

import re
from config.data_models.output_types import GromacsOutputs, GromacsPaths
from config.mdp_workflow_config import minim_workflow, polymer_workflow
from modules.workflows.base_workflow import BaseWorkflow
from config.paths import TEMP_DIR
from modules.utils.shared.file_utils import delete_directory
from typing import Any, Dict, Tuple
import math
import csv
import logging
import os

from modules.cache_store.equilibriated_atomistic_polymer_cache import (
    EquilibriatedAtomisticPolymerCache,
)


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
equiibriated_atomistic_polymer_cache = EquilibriatedAtomisticPolymerCache()


class NoJoinedAtomisticPolymerWorkflow(BaseWorkflow):
    analysis_folder = "analysis"
    absolute_min_box_size: float = 5.0
    average_bond_length: float = 0.15
    safety_factor: float = 1.5
    standard_box_width: float = 7.0

    def __init__(
        self,
        input_dir: str,
        solvent_name: str,
        solvent_smiles: str,
        solvent_density: float,
        temperatures: List[float],
        output_dir: str,
        solvent_compressibility: float,
        csv_file_path: str,
        minim_workflow: FullEquilibrationWorkflow = minim_workflow,
        full_workflow: FullEquilibrationWorkflow = polymer_workflow,
        cache: EquilibriatedAtomisticPolymerCache = equiibriated_atomistic_polymer_cache,
        pos_ion_name: str = "NA",
        neg_ion_name: str = "CL",
        polymer_name: str = "UNL",
        sol_resname: str = "SOL",
        verbose: bool = True,
        cleanup: bool = True,
        confirm_temp_deletion: bool = True,
        box_incriments: float = 5,
    ):
        self.csv_file_path = f"{csv_file_path}.csv"
        self.solvent_smiles = solvent_smiles
        self.cleanup = cleanup
        self.confirm_temp_deletion = confirm_temp_deletion
        self.solvent = SolventGenerator(
            solvent_name=solvent_name,
            solvent_smiles=solvent_smiles,
            solvent_compressibility=solvent_compressibility,
            solvent_density=solvent_density,
            solvent_resname=sol_resname,
        ).generate_solvent()

        self.temperatures = temperatures
        self.directory_contents: List[Tuple[str, int, GromacsPaths]] = (
            self.parse_directory(input_dir)
        )
        self.output_dir = output_dir
        self.minim_workflow = minim_workflow
        self.full_workflow = full_workflow
        self.cache = cache
        self.pos_ion_name = pos_ion_name
        self.neg_ion_name = neg_ion_name
        self.polymer_name = polymer_name
        self.verbose = verbose

    @staticmethod
    def parse_directory(base_dir: str) -> List[Tuple[str, int, GromacsPaths]]:

        results = []
        pattern = re.compile(
            r"(.+)_(\d+)$"
        )  # Match monomer_smiles and n_units at the end

        for subdir in os.listdir(base_dir):

            subdir_path = os.path.join(base_dir, subdir)
            if not os.path.isdir(subdir_path):
                continue  # Skip files, only process directories

            match = pattern.match(subdir)
            if match:
                monomer_smiles, n_units = match.groups()
                n_units = int(n_units)  # Convert to integer

                gromacs_paths = GromacsPaths(
                    itp_path=(
                        os.path.join(subdir_path, "POLY_GMX.itp")
                        if os.path.exists(os.path.join(subdir_path, "POLY_GMX.itp"))
                        else None
                    ),
                    gro_path=(
                        os.path.join(subdir_path, "POLY_GMX.gro")
                        if os.path.exists(os.path.join(subdir_path, "POLY_GMX.gro"))
                        else None
                    ),
                    top_path=(
                        os.path.join(subdir_path, "POLY_GMX.top")
                        if os.path.exists(os.path.join(subdir_path, "POLY_GMX.top"))
                        else None
                    ),
                    posre_path=(
                        os.path.join(subdir_path, "POLY_GMX_posre.itp")
                        if os.path.exists(
                            os.path.join(subdir_path, "POLY_GMX_posre.itp")
                        )
                        else None
                    ),
                )

                results.append((monomer_smiles, n_units, gromacs_paths))

        return results

    def _run_simulations(
        self,
        temperature: float,
        num_units: int,
        monomer_smiles: List[str],
        parameterised_polymer: GromacsPaths,
    ) -> GromacsOutputs:
        box_size_nm = self._get_min_box_size(num_units=num_units)
        polymer_workflow = NoParametiserPolymerEquilibriationWorkflow(
            parametised_polymer=parameterised_polymer,
            monomer_smiles=monomer_smiles,
            num_units=num_units,
            solvent=self.solvent,
            box_size_nm=box_size_nm,
            temperatures=[temperature],
            output_dir=self.output_dir,
            minim_workflow=self.minim_workflow,
            full_workflow=self.full_workflow,
            cache=self.cache,
            pos_ion_name=self.pos_ion_name,
            neg_ion_name=self.neg_ion_name,
            polymer_name=self.polymer_name,
            verbose=self.verbose,
            cleanup_log=True,
            cleanup_temp=False,
            confirm_temp_deletion=False,
        )
        outputs = polymer_workflow.run()[0]
        n_values = polymer_workflow.actual_num_units
        final_output_dir = polymer_workflow.final_output_dir
        return outputs, n_values, final_output_dir

    def _get_min_box_size(self, num_units: int, box_incriments: float = 5):
        if num_units <= 30:
            box_width = self.standard_box_width
        if num_units > 30:
            end_to_end_length = num_units * self.average_bond_length
            safe_size = end_to_end_length * self.safety_factor
            box_width = math.ceil(safe_size / box_incriments) * box_incriments

        return [box_width, box_width, box_width]

    def _analyse_outputs(
        self, outputs: GromacsOutputs, output_dir: str, temperature: float
    ):
        analysis_dir = os.path.join(output_dir, self.analysis_folder)
        temperature_dir = os.path.join(analysis_dir, f"T_{temperature}")
        os.makedirs(temperature_dir, exist_ok=True)
        analyser = GromacsAnalyser(
            outputs=outputs,
            poly_resname=self.polymer_name,
            ion_resnames=[self.pos_ion_name, self.neg_ion_name],
            output_dir=temperature_dir,
        )
        Rg_mean, Rg_std = analyser.extract_radius_of_gyration()
        D = analyser.extract_diffusion_coefficient()
        SASA_mean, SASA_std = analyser.extract_sasa()
        E2E_mean, E2E_std = analyser.extract_end_to_end_distance()
        return Rg_mean, Rg_std, D, SASA_mean, SASA_std, E2E_mean, E2E_std

    def _run_per_temp(
        self,
        temperature: float,
        num_units: int,
        monomer_smiles: List[str],
        parameterised_polymer: GromacsPaths,
    ) -> str:
        outputs, n_units, final_output_dir = self._run_simulations(
            temperature=temperature,
            num_units=num_units,
            monomer_smiles=monomer_smiles,
            parameterised_polymer=parameterised_polymer,
        )
        Rg_mean, Rg_std, D, SASA_mean, SASA_std, E2E_mean, E2E_std = (
            self._analyse_outputs(
                outputs=outputs, output_dir=final_output_dir, temperature=temperature
            )
        )
        row_data = {
            "solvent_name": self.solvent.name,
            "solvent_smiles": self.solvent_smiles,
            "compressibility": self.solvent.compressibility,
            "monomer_smiles_list": ";".join(self.monomer_smiles),
            "N": n_units,
            "T": temperature,
            "Rg_mean": Rg_mean,
            "Rg_std": Rg_std,
            "Diffusion_Coefficient": D,
            "SASA_mean": SASA_mean,
            "SASA_std": SASA_std,
            "E2E_mean": E2E_mean,
            "E2E_std": E2E_std,
        }
        logger.info(f"Row data: {row_data}")
        csv = self._write_csv_row(row_data)
        return csv

    def run(self):
        for monomer_smiles, num_units, parameterised_polymer in self.directory_contents:
            for temperature in self.temperatures:
                self._run_per_temp(
                    temperature=temperature,
                    num_units=num_units,
                    monomer_smiles=monomer_smiles,
                    parameterised_polymer=parameterised_polymer,
                )
            if self.cleanup_temp:
                delete_directory(
                    TEMP_DIR, verbose=self.verbose, confirm=self.confirm_temp_deletion
                )

        return self.csv_file_path

    def _write_csv_row(self, row_data: Dict[str, Any]) -> str:
        with open(self.csv_file_path, "a") as f:
            writer = csv.DictWriter(f, fieldnames=row_data.keys())
            writer.writerow(row_data)

        logger.info(f"Row appended to {self.csv_file_path}")
        return self.csv_file_path

    def clear_past_csv(self, confirm: bool = True):
        if self.csv_file_path.exists():
            if confirm:
                confirm = input(
                    f"Are you sure you want to delete {self.csv_file_path}? (y/n) "
                )
            if confirm == "y":
                os.remove(self.csv_file_path)
                logger.info(f"Deleted {self.csv_file_path}")
