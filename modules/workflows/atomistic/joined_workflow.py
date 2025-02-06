from modules.rdkit.solvent_generator import SolventGenerator
from modules.workflows.atomistic.polymer_equilibriator import (
    PolymerEquilibriationWorkflow,
)
from modules.gromacs.analyser import GromacsAnalyser
from typing import List, Any, Dict
from modules.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)

from config.data_models.output_types import GromacsOutputs
from config.mdp_workflow_config import minim_workflow, polymer_workflow
from modules.workflows.base_workflow import BaseWorkflow
from config.paths import TEMP_DIR
from modules.utils.shared.file_utils import delete_directory
from typing import Any, Dict
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


class JoinedAtomisticPolymerWorkflow(BaseWorkflow):
    analysis_folder = "analysis"
    absolute_min_box_size: float = 5.0
    average_bond_length: float = 0.15
    safety_factor: float = 1.5
    standard_box_width: float = 7.0

    def __init__(
        self,
        monomer_smiles: List[str],
        num_units: int,
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
<<<<<<< HEAD
        polymer_mol_name: str = "POLY",
=======
>>>>>>> 91758eb (cleaned up)
        sol_resname: str = "SOL",
        verbose: bool = True,
        cleanup: bool = True,
        confirm_temp_deletion: bool = True,
        box_incriments: float = 5,
    ):
<<<<<<< HEAD
        self.polymer_mol_name = polymer_mol_name
        self.cleanup_temp = cleanup
=======
>>>>>>> 91758eb (cleaned up)
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
        self.num_units = num_units
        self.monomer_smiles = monomer_smiles
        self.output_dir = output_dir
        self.minim_workflow = minim_workflow
        self.full_workflow = full_workflow
        self.cache = cache
        self.pos_ion_name = pos_ion_name
        self.neg_ion_name = neg_ion_name
        self.polymer_name = polymer_name
        self.verbose = verbose
<<<<<<< HEAD
        self.data = None
=======
>>>>>>> 91758eb (cleaned up)
        self.box_size_nm = self._get_min_box_size(box_incriments=box_incriments)

    def _run_simulations(self, temperature: float) -> GromacsOutputs:
        polymer_workflow = PolymerEquilibriationWorkflow(
            monomer_smiles=self.monomer_smiles,
            num_units=self.num_units,
            solvent=self.solvent,
            box_size_nm=self.box_size_nm,
            temperatures=[temperature],
            output_dir=self.output_dir,
            minim_workflow=self.minim_workflow,
            full_workflow=self.full_workflow,
            cache=self.cache,
            pos_ion_name=self.pos_ion_name,
            neg_ion_name=self.neg_ion_name,
<<<<<<< HEAD
            polymer_name=self.polymer_mol_name,
=======
            polymer_name=self.polymer_name,
>>>>>>> 91758eb (cleaned up)
            verbose=self.verbose,
            cleanup_log=True,
            cleanup_temp=False,
            confirm_temp_deletion=False,
        )
        outputs = polymer_workflow.run()[0]
        n_values = polymer_workflow.actual_num_units
        final_output_dir = polymer_workflow.final_output_dir
        return outputs, n_values, final_output_dir

    def _get_min_box_size(self, box_incriments: float = 5):
        if self.num_units <= 30:
            box_width = self.standard_box_width
        if self.num_units > 30:
            end_to_end_length = self.num_units * self.average_bond_length
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

    def _run_per_temp(self, temperature: float) -> str:
        outputs, n_units, final_output_dir = self._run_simulations(
            temperature=temperature
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
<<<<<<< HEAD
        self.data = row_data
=======
>>>>>>> 91758eb (cleaned up)
        csv = self._write_csv_row(row_data)
        return csv

    def run(self):
        for temperature in self.temperatures:
            self._run_per_temp(temperature=temperature)
        if self.cleanup_temp:
            delete_directory(
<<<<<<< HEAD
                TEMP_DIR, verbose=self.verbose, confirm=False
=======
                TEMP_DIR, verbose=self.verbose, confirm=self.confirm_temp_deletion
>>>>>>> 91758eb (cleaned up)
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
