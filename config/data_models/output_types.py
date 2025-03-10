from dataclasses import dataclass
from typing import Optional, List
from modules.cg_mappers.martini_index_generator import MARTINIIndexGenerator
from modules.cg_mappers.martini_map_generator import MARTINIMapGenerator
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator

from dataclasses import dataclass


@dataclass
class GromacsPaths:
    """
    A class to store the paths of the acpype files.
    """

    itp_path: Optional[str] = None
    gro_path: Optional[str] = None
    top_path: Optional[str] = None
    posre_path: Optional[str] = None

    def to_list(self) -> List[Optional[str]]:
        """
        Converts the GromacsPaths instance into a list of file paths.

        :return: A list of file paths.
        :rtype: List[Optional[str]]
        """
        print("!!!!!!!!!!!!!!!4")
        return [self.itp_path, self.gro_path, self.top_path, self.posre_path]


@dataclass
class GromacsOutputs:
    itp: Optional[str] = None
    gro: Optional[str] = None
    top: Optional[str] = None
    edr: Optional[str] = None
    trr: Optional[str] = None
    log: Optional[str] = None
    tpr: Optional[str] = None
    xtc: Optional[str] = None


class MARTINIMaps:
    """
    Manages the generation of MARTINI .map and .ndx files.
    """

    def __init__(self, polymer: BasePolymerGenerator, file_name: str, output_dir: str):
        self.polymer = polymer
        self.file_name = file_name
        self.output_dir = output_dir

    @property
    def map_file(self) -> str:
        """
        Generates the MARTINI .map file only when accessed.
        """
        return MARTINIMapGenerator(self.polymer).create_map(
            self.file_name, self.output_dir
        )

    @property
    def ndx_file(self) -> str:
        """
        Generates the MARTINI .ndx file only when accessed.
        """
        return MARTINIIndexGenerator(self.polymer).create_map(
            self.file_name, self.output_dir
        )
