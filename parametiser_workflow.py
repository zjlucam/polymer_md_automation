from modules.workflows.separated.parametiser.solvent_csv import SolventCSVParametiser
from modules.workflows.separated.parametiser.polymer_list import PolymerListParametiser

<<<<<<< HEAD
#solvent_generator = SolventCSVParametiser(
#    solvent_csv_path="input_data/solvent_data.csv", output_dir="test_solvent"
#).run()
=======
solvent_generator = SolventCSVParametiser(
    solvent_csv_path="input_data/solvent_data.csv", output_dir="test_solvent"
).run
>>>>>>> 91758eb (cleaned up)

from input_data.monomer_smiles import monomer_smiles_list

PolymerListParametiser = PolymerListParametiser(
    full_smiles_list=monomer_smiles_list,
    output_dir="test_polymer",
<<<<<<< HEAD
    num_units_list=[20],
=======
    num_units_list=[],
>>>>>>> 91758eb (cleaned up)
).run()
