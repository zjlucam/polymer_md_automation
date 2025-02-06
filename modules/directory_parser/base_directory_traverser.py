from abc import ABC, abstractmethod
from typing import List
import os


class BaseDirectoryParser(ABC):

    def __init__(self, base_dir: str):
        self.base_dir = base_dir

    def parse_directory(self) -> List:

        results = []
        for subdir in os.listdir(self.base_dir):
            subdir_path = os.path.join(self.base_dir, subdir)
            if os.path.isdir(subdir_path):
                result = self.parse_subdirectory(subdir_path)
                if result:
                    results.append(result)
        return results

    @abstractmethod
    def parse_subdirectory(self, subdir_path: str):

        pass
