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

    manager.run()
