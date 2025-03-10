import pickle
import os
import logging
from typing import Any
from config.paths import MAIN_CACHE_DIR
from modules.cache_store.base_cache import BaseCache

logger = logging.getLogger(__name__)


class PickleCache(BaseCache):
    """
    A cache for storing and retrieving objects using Pickle.
    """

    def __init__(self, name: str, cache_dir=MAIN_CACHE_DIR):
        """
        Initializes the PickleCache.

        :param cache_dir: Directory to store cached objects.
        """
        super().__init__(cache_dir=cache_dir, cache_name=f"picklecache_{name}")

    def _serialize(self, data: Any) -> str:
        """
        Serializes data using Pickle and stores the file path in the cache.

        :param data: The data to serialize.
        :return: The file path where the serialized data is stored.
        """
        hash_key = self._generate_hash(data)
        file_path = os.path.join(self.cache_dir, f"{hash_key}.pkl")

        with open(file_path, "wb") as file:
            pickle.dump(data, file)

        return file_path

    def _deserialize(self, file_path: str) -> Any:
        """
        Deserializes Pickle data from a file.

        :param file_path: Path to the Pickle file.
        :return: The deserialized object.
        """
        if not os.path.exists(file_path):
            return None

        with open(file_path, "rb") as file:
            return pickle.load(file)

    def _generate_hash(self, data: Any) -> str:
        """
        Generates a hash key based on the object's Pickle serialization.

        :param data: The data object to hash.
        :return: A unique hash string.
        """
        return str(hash(pickle.dumps(data)))

    def store_object(self, key: str, data: Any):
        """
        Stores an object in the cache.

        :param key: Unique key for retrieving the object later.
        :param data: The object to store.
        """
        file_path = self._serialize(data)
        self.store(key, file_path)
        logger.info(f"Stored object in cache at {file_path}")

    def retrieve_object(self, key: str) -> Any:
        """
        Retrieves an object from the cache.

        :param key: The key for the object.
        :return: The cached object if available, otherwise None.
        """
        if not self.has_key(key):
            return None

        file_path = self.retrieve(key)
        return self._deserialize(file_path)
