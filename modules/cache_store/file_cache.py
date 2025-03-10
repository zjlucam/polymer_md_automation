from modules.cache_store.base_cache import BaseCache
from pathlib import Path
from config.paths import MAIN_CACHE_DIR
from typing import Any, Optional
import os
import logging


logger = logging.getLogger(__name__)


class FileCache(BaseCache):
    """Caches file paths for quick lookup."""

    def __init__(self, name: str, cache_dir=MAIN_CACHE_DIR):
        """
        Initializes the PickleCache.

        :param cache_dir: Directory to store cached objects.
        """
        super().__init__(
            cache_name=f"file_cache_{name}",
            cache_dir=cache_dir,
        )

    def _serialize(self, data: Any) -> str:
        """
        Serializes the file path as a string.

        :param data: The file path to serialize.
        :return: The serialized file path.
        """
        return str(Path(data).resolve())  # Ensure absolute path as a string

    def _deserialize(self, data: Any) -> Optional[str]:
        """
        Deserializes the stored file path.

        :param data: The stored data.
        :return: The file path as a string, or None if invalid.
        """
        return str(Path(data).resolve()) if data else None

    def store_object(self, key: str, file_path: str):
        """
        Stores a file path in the cache.

        :param key: The key associated with the file path.
        :param file_path: The file path to store.
        """
        if not os.path.exists(file_path):
            raise ValueError(f"File does not exist: {file_path}")

        file_path = self._serialize(file_path)  # Convert to absolute path string
        self.store(key, file_path)  # Use BaseCache store method
        logger.info(f"Stored file path in cache: {file_path}")

    def retrieve_object(self, key: str) -> Optional[str]:
        """
        Retrieves a file path from the cache.

        :param key: The key associated with the file path.
        :return: The cached file path if found, otherwise None.
        """
        file_path = self.retrieve(key)  # Use BaseCache retrieve method
        return self._deserialize(file_path) if file_path else None

    def file_exists(self, key: str) -> bool:
        """
        Checks if a file exists at the cached path.

        :param key: The key associated with the file path.
        :return: True if the file exists, False otherwise.
        """
        file_path = self.retrieve_object(key)
        return file_path is not None and os.path.exists(file_path)
