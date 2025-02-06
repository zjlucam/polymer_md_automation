<<<<<<< HEAD
from input_data.monomer_smiles import monomer_smiles_list
from simulation_manager import PolymerSimulationManager

if __name__ == "__main__":
    manager = PolymerSimulationManager(
        solvent_csv="input_data/solvent_data.csv",
        monomer_smiles=monomer_smiles_list,
        num_units=[20, 25, 15],
        temperatures=[280, 298, 348],
        output_dir="outputs_test_run",
        csv_file_path="output_2_4.csv",
    )

=======
from simulation_manager import SimulationManager
import datetime
import os


job_id = os.getenv("SLURM_JOB_ID", "local_run")
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
timestamp = datetime.datetime.now().strftime("%m-%d_%H")

if __name__ == "__main__":
    manager = SimulationManager(
        parameterised_poly_dir="test_sol",
        parameterised_sol_dir="test_poly",
        temperatures=[298],
        output_dir=f"output_{job_id}_{timestamp}",
        identifying_tag=job_id,
        polymer_start_idx=5,
        solvent_end_idx=5,
        polymer_end_idx=25,
        solvent_end_idx=10,
    )
>>>>>>> 91758eb (cleaned up)
    manager.run()
