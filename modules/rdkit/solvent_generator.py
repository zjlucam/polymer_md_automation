from modules.rdkit.base_molecule_generator import BaseMoleculeGenerator
from config.paths import SOLVENT_PDB_DIR
from modules.cache_store.file_cache import FileCache
from config.data_models.solvent import Solvent
import re
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

solvent_pdb_cache = FileCache(name="solvent_pdb_cache", cache_dir=SOLVENT_PDB_DIR)


class SolventGenerator(BaseMoleculeGenerator):
    def __init__(
        self,
        solvent_name: str,
        solvent_smiles: str,
        solvent_compressibility: str,
        solvent_density: float,
        solvent_resname: str = "SOL",
    ):
        self.sol_name = solvent_name
        self.sol_density = solvent_density
        self.sol_smiles = solvent_smiles
        self.sol_compressibility = solvent_compressibility
        self.sol_molecule = self._convert_to_rdkit_molecule(self.sol_smiles)
        self.sol_molar_mass = self._get_molar_mass(self.sol_molecule)
        self.cache = solvent_pdb_cache
        self.cache_key = self._generate_filename()
        self.sol_resname = solvent_resname

    def _generate_solvent_mol(self):
        mol = self._convert_to_rdkit_molecule(self.sol_smiles)
        mol = self._finalise_molecule(mol)
        return mol

    def _rename_res_in_pdb(self, pdb_path: str) -> str:
        with open(pdb_path, "r") as file:
            pdb_lines = file.readlines()

        with open(pdb_path, "w") as file:
            for line in pdb_lines:
                if " UNL " in line:
                    line = line.replace(" UNL ", f" {self.sol_resname.ljust(3)} ")
                file.write(line)
        return pdb_path

    @staticmethod
    def sanitize_filename(smiles: str) -> str:
        """
        Converts a SMILES string into a valid filename by replacing
        special characters with underscores.

        :param smiles: SMILES string to be sanitized
        :return: Sanitized filename string
        """
        return re.sub(r"[^\w\-_]", "_", smiles)

    def _generate_filename(self):
        return self.sanitize_filename(self.sol_smiles)

    def check_solvent_pdb_cache(self, cache_key: str):
        if self.cache.has_key(cache_key):
            logger.info(f"Retrieving solvent pdb from cache: {cache_key}")
            return self.cache.retrieve_object(cache_key)
        return None

    def _get_solvent_pdb(self, output_dir: str = SOLVENT_PDB_DIR) -> str:
        cache_key = self.cache_key
        cached_mol_pdb = self.check_solvent_pdb_cache(cache_key=cache_key)
        if cached_mol_pdb:
            logger.info(f"Retrieved solvent pdb from cache: {cache_key}")
            return cached_mol_pdb
        logger.info("No cached solvent pdb found, generating")
        logger.info(f"Generating solvent pdb for {self.sol_name}")
        mol = self._generate_solvent_mol()
        pdb_path = self._save_as_pdb(mol, output_dir, self._generate_filename())
        self.cache.store_object(cache_key, pdb_path)
        logger.info(f"Stored solvent pdb in cache: {cache_key}")
        return pdb_path

    def generate_solvent(self, output_dir: str = SOLVENT_PDB_DIR) -> Solvent:
        pdb_path = self._get_solvent_pdb(output_dir=output_dir)
        pdb_path = self._rename_res_in_pdb(pdb_path)
        solvent = Solvent(
            name=self.sol_name,
            density=self.sol_density,
            molecular_weight=self.sol_molar_mass,
            pdb_path=pdb_path,
            compressibility=self.sol_compressibility,
        )

        return solvent
