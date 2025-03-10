from modules.workflows.separated.parametiser.solvent import SolventParametiser
import pandas as pd
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

class SolventCSVParametiser:

    def __init__(
        self,
        solvent_csv_path: str,
        output_dir: str,
        solvent_resname: str = "SOL",
        verbose: bool = False,
    ):
        self.output_dir = output_dir
        self.verbose = verbose
        self.solvent_resname = solvent_resname
        self.solvent_df = pd.read_csv(solvent_csv_path)
        self.solvent_csv_path = solvent_csv_path

    def run(self):
        for _, solvent in self.solvent_df.iterrows():
            try:
                solvent_name = solvent["name"]
                solvent_smiles = solvent["SMILES"]
                solvent_compressibility = solvent["compressibility"]
                solvent_density = solvent["density"]

                generator = SolventParametiser(
                    solvent_name=solvent_name,
                    solvent_smiles=solvent_smiles,
                    solvent_compressibility=solvent_compressibility,
                    solvent_density=solvent_density,
                    output_dir=self.output_dir,
                    verbose=self.verbose,
                    solvent_resname=self.solvent_resname,
                )
                generator.run()
            except ValueError as e:
                logger.error(f"ValueError in {solvent_name}: {e}")