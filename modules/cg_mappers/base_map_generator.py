import xml.etree.ElementTree as ET
from typing import List, Dict
import xml.etree.ElementTree as ET
from typing import List, Tuple, Dict, Optional
import xml.dom.minidom
import pandas as pd
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.data_handler import DataHandler
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from modules.utils.shared.file_utils import check_directory_exists
from abc import ABC, abstractmethod
import logging
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BaseMapGenerator(ABC):
    """
    Generates a VOTCA-compatible XML mapping file for coarse-grained simulations.

    Attributes:
        bead_mappings (List[Dict]): Coarse-grained bead mappings.
        bonds (List[Tuple[int, int]]): List of coarse-grained bonds.
        angles (List[Tuple[int, int, int]]): List of coarse-grained angles.
    """

    map_extension: str

    def __init__(self, polymer: BasePolymerGenerator):
        """
        Initializes the BaseMapGenerator with bead mappings, bonds, and angles.

        Args:
            polymer (BasePolymerGenerator): The polymer object containing CG mapping information.
        """
        super().__init__()
        self.bead_mappings = polymer.cg_map
        self.bonds = polymer.cg_bonds if polymer.cg_bonds else []
        self.angles = polymer.cg_angles if polymer.cg_angles else []

    def _generate_filepath(self, filename: str, output_dir: Optional[str]) -> str:
        """
        Generates a filepath for the mapping file.

        Args:
            filename (str): The name of the file.
            output_dir (Optional[str]): The directory to save the file.

        Returns:
            str: The full file path.
        """
        if output_dir:
            check_directory_exists(output_dir, make_dirs=True)
            return f"{output_dir}/{filename}.{self.map_extension}"
        return f"{filename}.{self.map_extension}"

    @abstractmethod
    def _generate_mapping(self, start_index: Optional[int] = None) -> str:
        """
        Abstract method to generate the mapping content.

        Args:
            start_index (Optional[int]): The starting index for atom numbering.

        Returns:
            str: The generated mapping content.
        """
        pass

    def _save_mapping(
        self,
        content: str,
        filename: str,
        output_dir: Optional[str],
        start_index: Optional[int] = None,
    ) -> str:
        """
        Saves the generated mapping content to a file.

        Args:
            content (str): The content to save.
            filename (str): The name of the file.
            output_dir (Optional[str]): The directory to save the file.
            start_index (Optional[int]): The starting index for atom numbering.

        Returns:
            str: The saved file path.
        """
        file_path = self._generate_filepath(filename, output_dir)
        with open(file_path, "w") as f:
            f.write(content)
        return file_path

    def create_map(
        self,
        filename: str,
        output_dir: Optional[str] = None,
        start_index: Optional[int] = None,
    ) -> str:
        """
        Creates the mapping file by generating and saving the mapping content.

        Args:
            filename (str): The name of the mapping file.
            output_dir (Optional[str]): The directory to save the file.
            start_index (Optional[int]): The starting index for atom numbering.

        Returns:
            str: The saved file path.
        """
        content = self._generate_mapping(start_index=start_index)
        return self._save_mapping(content, filename, output_dir, start_index)

    @staticmethod
    def reformat_atom_name(atom_name: str, start_index: Optional[int] = None) -> str:
        """
        Reformats atom names based on the start_index.
        - If `start_index` is None, preserves original naming (no C0, H0, etc.).
        - If `start_index` is provided, shifts the numbering accordingly.

        Args:
            atom_name (str): The original atom name.
            start_index (Optional[int]): The index offset.

        Returns:
            str: The reformatted atom name.
        """
        match = re.match(r"^([A-Za-z]+)(\d*)$", atom_name)
        if not match:
            raise ValueError(f"Unexpected atom name format: {atom_name}")

        atom_type, index = match.groups()

        if index == "" and start_index is None:
            return atom_type  # Preserve original format (C -> C, H -> H)

        index = int(index) if index else 0  # Convert to integer, default to 0 if empty
        new_index = index + (start_index or 0)

        if start_index is None and new_index == 0:
            return atom_type  # Avoid adding "0" when start_index is None

        return f"{atom_type}{new_index}"
