<<<<<<< HEAD
import itertools
import logging
import os
import pandas as pd
from typing import List, Tuple
from modules.workflows.atomistic.joined_workflow import JoinedAtomisticPolymerWorkflow
import random
import multiprocessing
import json


# Set up logging
=======
from modules.directory_parser.polymer_directory_parser import PolymerDirectoryParser
from modules.directory_parser.solvent_directory_parser import SolventDirectoryParser
from modules.workflows.separated.gromacs.joined import JoinedAtomisticPolymerWorkflow
from modules.utils.shared.file_utils import check_directory_exists
import logging
import os
import pandas as pd
from typing import List

>>>>>>> 91758eb (cleaned up)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


<<<<<<< HEAD
class PolymerSimulationManager:
    error_file_path = "error.csv"
    random_seed = 42
    def __init__(
        self,
        solvent_csv: str,
        monomer_smiles: List[str],
        output_dir: str,
        progress_file: str = "progress",
        csv_file_path: str = "outputs",
        num_units: List[int] = [5, 10, 20],
        temperatures: List[int] = [280, 298, 346],
    ):
        """
        Initializes the simulation manager.
        :param solvent_csv: Path to the CSV file containing solvent properties.
        :param monomer_smiles: List of monomer SMILES strings.
        :param output_dir: Directory where results should be stored.
        :param progress_file: File for tracking progress.
        """
        self.solvent_df = pd.read_csv(solvent_csv)
        self.monomer_smiles = monomer_smiles
        self.num_units = num_units
        self.temperatures = temperatures
        self.output_dir = output_dir
        self.progress_file = f"{progress_file}.csv"
        self.csv_file_path = csv_file_path
        self._initialise_progress_file()
        self.completed_jobs = self._load_progress()

    def _initialise_progress_file(self):
        if not os.path.exists(self.progress_file):
            pd.DataFrame(columns =["job_id"]).to_csv(self.progress_file, index = False)
            logger.info(f"Created progress file at {self.progress_file}")

    def _load_progress(self):
        """
        Load completed jobs from progress file to avoid re-running.
        """
        if os.path.exists(self.progress_file):
            try:
                progress_df = pd.read_csv(self.progress_file)
                return set(progress_df["job_id"].astype(str))
            except Exception as e:
                logger.error(f"Failed to read progress file {self.progress_file}:{e}")
            return set()
        return set()

    def _save_progress(self, job_id: str):
        """
        Save the completed job to progress file.
        """
        progress_entry = pd.DataFrame([[job_id]], columns=["job_id"])
        progress_entry.to_csv(self.progress_file, mode = "a", header = not os.path.exists(self.progress_file), index = False)
        logger.info(f"Progress saved for {job_id}")

    def _log_error(self, job_id: str, error_message: str):
        error_entry = pd.DataFrame([[job_id, error_message]], columns=["Job ID", "Error"])
        error_entry.to_csv(self.error_file_path, mode = "a", header = False, index = False)
        

    def _generate_combinations(self):
        """
        Generate all single-monomer, two-monomer, and three-monomer combinations with solvents.
        """
        
        all_combinations = []
        for num_units in self.num_units:
            for temp in self.temperatures:
                for _, solvent in self.solvent_df.iterrows():
                    solvent_name = solvent["name"]
                    solvent_smiles = solvent["SMILES"]
                    solvent_compressibility = float(solvent["compressibility"])
                    solvent_density = float(solvent["density"])

                    # Single-monomer cases
                    for monomer in self.monomer_smiles:
                        all_combinations.append(
                            (
                                [monomer],
                                solvent_name,
                                solvent_smiles,
                                solvent_density,
                                solvent_compressibility,
                                temp,
                                num_units
                            )
                        )

                    # Two-monomer copolymers (order-independent, no duplicates)
                    for monomer_pair in itertools.combinations(self.monomer_smiles, 2):
                        all_combinations.append(
                            (
                                list(monomer_pair),
                                solvent_name,
                                solvent_smiles,
                                solvent_density,
                                solvent_compressibility,
                                temp,
                                num_units
                            )
                        )

                    # Three-monomer copolymers
                    for monomer_triplet in itertools.combinations(self.monomer_smiles, 3):
                        all_combinations.append(
                            (
                                list(monomer_triplet),
                                solvent_name,
                                solvent_smiles,
                                solvent_density,
                                solvent_compressibility,
                                temp,
                                num_units
                            )
                        )

        return all_combinations
    
    def run(self, timeout=3200):
        """
        Iterates through all polymer-solvent-temperature combinations and executes the workflow.
        Includes a timeout mechanism to avoid stuck simulations.
        """
        self.completed_jobs = self._load_progress()
        combinations = self._generate_combinations()
        random.seed(self.random_seed)
        random.shuffle(combinations)

        for (
            monomer_list,
            solvent_name,
            solvent_smiles,
            solvent_density,
            solvent_compressibility,
            temp,
            num_units
        ) in combinations:
            job_id = f"{'_'.join(monomer_list)}_{solvent_name}_{num_units}_{temp}"

            if job_id in self.completed_jobs:
                logger.info(f"Skipping {job_id}, already completed.")
                continue

            try:
                logger.info(f"Running simulation for: {job_id}")

                # Run workflow in a separate process with a timeout
                process = multiprocessing.Process(
                    target=self._run_workflow, 
                    args=(monomer_list, solvent_name, solvent_smiles, 
                          solvent_density, solvent_compressibility, temp, num_units, job_id)
                )
                process.start()
                process.join(timeout)  # Wait for workflow to complete or timeout
                
                if process.is_alive():
                    logger.warning(f"Timeout exceeded for {job_id}, terminating process.")
                    process.terminate()
                    process.join()
                    self._log_error(job_id=job_id, error_message="Simulation timed out.")

                else:
                    # Log progress if completed successfully
                    self._save_progress(job_id)

            except ValueError as e:
                error_message = f"ValueError in {job_id}: {e}"
                self._log_error(job_id=job_id, error_message=error_message)
            except Exception as e:
                error_message = f"Unexpected error in {job_id}: {e}"
                self._log_error(job_id=job_id, error_message=error_message)

        logger.info("All simulations completed!")

    def _run_workflow(self, monomer_list, solvent_name, solvent_smiles, 
                      solvent_density, solvent_compressibility, temp, num_units, job_id):
        """
        Wrapper function to run workflow. This runs in a separate process.
        """
        try:
            workflow = JoinedAtomisticPolymerWorkflow(
                monomer_smiles=monomer_list,
                num_units=num_units,
                solvent_name=solvent_name,
                solvent_smiles=solvent_smiles,
                solvent_density=solvent_density,
                temperatures=[temp],
                output_dir=self.output_dir,
                solvent_compressibility=solvent_compressibility,
                csv_file_path=self.csv_file_path,
            )

            workflow.run()  # Runs the simulation
            
            data = workflow.data
            
                    # Prepare data for CSV
            output_file = "output_backup.csv"
            row = {
                "monomer_list": str(monomer_list),
                "solvent_name": solvent_name,
                "solvent_smiles": solvent_smiles,
                "solvent_density": solvent_density,
                "solvent_compressibility": solvent_compressibility,
                "temperature": temp,
                "num_units": num_units,
                "job_id": job_id,
                "data": json.dumps(data)  # Convert data to string if it's not scalar
            }

            # Convert to DataFrame
            df = pd.DataFrame([row])

            # Append or create CSV
            if not os.path.exists(output_file):
                df.to_csv(output_file, index=False)
            else:
                df.to_csv(output_file, mode='a', header=False, index=False)

        except Exception as e:
            logger.error(f"Error in {job_id}: {e}")
=======
class SimulationManager:
    error_csv_filename = "error_log"
    output_csv_filename = "outputs"  # .csv

    def __init__(
        self,
        parameterised_poly_dir: str,
        parameterised_sol_dir: str,
        temperatures: List[float],
        output_dir: str,
        identifying_tag: str = "",
        solvent_start_idx: int = 0,
        solvent_end_idx: int = None,
        polymer_start_idx: int = 0,
        polymer_end_idx: int = None,
    ):
        self.parameterised_solvents = SolventDirectoryParser(
            parameterised_sol_dir
        ).parse_directory()
        self.parameterised_polymers = PolymerDirectoryParser(
            parameterised_poly_dir
        ).parse_directory()

        # truncating so that I can split HPC submissions into multiple jobs
        solvent_length = len(self.parameterised_solvents)
        if solvent_end_idx is None or solvent_end_idx > solvent_length:
            solvent_end_idx = solvent_length

        polymer_length = len(self.parameterised_polymers)
        if polymer_end_idx is None or polymer_end_idx > polymer_length:
            polymer_end_idx = polymer_length

        self.parameterised_solvents = self.parameterised_solvents[
            solvent_start_idx:solvent_end_idx
        ]
        self.parameterised_polymers = self.parameterised_polymers[
            polymer_start_idx:polymer_end_idx
        ]
        self.output_dir = output_dir
        check_directory_exists(output_dir, make_dirs=True)
        self.temperatures = temperatures
        self.error_file_path = os.path.join(
            output_dir, f"{self.error_csv_filename}_{identifying_tag}.csv"
        )

    def _get_sanitised_monomer_smiles(self, monomer_smiles: List[str]) -> str:
        return "_".join(monomer_smiles)

    def run(self):
        for (
            parameterised_solvent,
            solvent,
            solvent_smiles,
        ) in self.parameterised_solvents:
            for (
                parameterised_polymer,
                monomer_smiles,
                n_units,
            ) in self.parameterised_polymers:
                for temp in self.temperatures:

                    job_id = f"{solvent.name}_{self._get_sanitised_monomer_smiles(monomer_smiles=monomer_smiles)}_{n_units}_{temp}"
                    try:
                        workflow = JoinedAtomisticPolymerWorkflow(
                            parameterised_polymer=parameterised_polymer,
                            monomer_smiles=monomer_smiles,
                            num_units=n_units,
                            solvent=solvent,
                            solvent_smiles=solvent_smiles,
                            parameterised_solvent=parameterised_solvent,
                            temperature=temp,
                            output_dir=self.output_dir,
                            csv_file_path=self.output_csv_filename,
                        ).run()
                        workflow.run()

                    except ValueError as e:
                        error_message = f"ValueError in {job_id}: {e}"
                        logger.error(error_message)
                        self._log_error(job_id, str(e))

                    except Exception as e:
                        error_message = f"Unexpected error in {job_id}: {e}"
                        logger.error(error_message)
                        self._log_error(job_id, str(e))

        logger.info("All simulations completed!")

    def _log_error(self, job_id: str, error_message: str):

        error_entry = pd.DataFrame(
            [[job_id, error_message]], columns=["Job ID", "Error"]
        )
        error_entry.to_csv(self.error_file_path, mode="a", header=False, index=False)
>>>>>>> 91758eb (cleaned up)
