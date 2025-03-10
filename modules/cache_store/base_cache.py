import os
import json
import pickle
import hashlib
from pathlib import Path
from abc import ABC, abstractmethod
from modules.utils.shared.file_utils import check_directory_exists
from config.paths import MAIN_CACHE_DIR
from typing import Any, Dict, Optional

import os
import json
import hashlib
import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


class BaseCache(ABC):
    """
    A base class for caching various types of objects.
    """

    def __init__(
        self,
        cache_name: str,
        cache_dir: str = MAIN_CACHE_DIR,
    ):
        check_directory_exists(directory_path=cache_dir, make_dirs=True)
        self.cache_dir = os.path.abspath(cache_dir)
        self.cache_index_path = os.path.join(self.cache_dir, f"{cache_name}_index.json")
        self.cache_index = self._load_cache_index()

    def _load_cache_index(self) -> Dict[str, Any]:
        if os.path.exists(self.cache_index_path):
            with open(self.cache_index_path, "r") as file:
                return json.load(file)
        return {}

    def _save_cache_index(self):
        """Saves the cache index to a JSON file."""
        with open(self.cache_index_path, "w") as file:
            json.dump(self.cache_index, file, indent=4)

    def has_key(self, key: str) -> bool:
        """
        Checks if a given key exists in the cache.

        :param key: The key to check.
        :return: True if the key exists, False otherwise.
        """
        return key in self.cache_index

    def store(self, key: str, data: Any):
        """
        Stores a value in the cache.

        :param key: The key under which to store the data.
        :param data: The data to store.
        """
        self.cache_index[key] = self._serialize(data)
        self._save_cache_index()

    def retrieve(self, key: str) -> Optional[Any]:
        """
        Retrieves a value from the cache.

        :param key: The key to retrieve.
        :return: The cached data if available, otherwise None.
        """
        return (
            self._deserialize(self.cache_index.get(key))
            if key in self.cache_index
            else None
        )

    def clear_cache(self):
        """Clears all cache data."""
        self.cache_index = {}
        self._save_cache_index()
