from modules.cache_store.pickle_cache import PickleCache
from config.paths import MAIN_CACHE_DIR
from config.data_models.solvent import Solvent
from typing import List


class EquilibriatedAtomisticPolymerCache(PickleCache):

    def __init__(self, cache_dir=MAIN_CACHE_DIR):
        super().__init__(name="polymer_cache", cache_dir=cache_dir)

    def get_cache_key(
        self,
        solvent: Solvent,
        monomer_smiles: List[str],
        num_units: float,
        temperature: float,
    ):
        monomer_smiles_str = "_".join(monomer_smiles)
        cache_key = f"{solvent.name}_{solvent.compressibility}_{monomer_smiles_str}_{num_units}_{temperature}"
        return cache_key
