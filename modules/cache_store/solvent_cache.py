from modules.cache_store.pickle_cache import PickleCache
from config.paths import MAIN_CACHE_DIR
from config.data_models.solvent import Solvent


class SolventCache(PickleCache):

    def __init__(self, cache_dir=MAIN_CACHE_DIR):
        super().__init__(name="solvent_cache", cache_dir=cache_dir)

    def get_cache_key(self, solvent: Solvent, temperature: float):
        return f"{solvent.name}_{solvent.compressibility}_{temperature}"
