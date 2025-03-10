from dataclasses import dataclass
from typing import Optional, List


@dataclass
class AcpypeOutputConfig:
    """
    A class to store the configuration for the acpype output.
    """

    itp: bool
    gro: bool
    top: bool
    posre: bool

    ITP_FILE_NAME_FORMAT: str = "{molecule_name}_GMX.itp"
    GRO_FILE_NAME_FORMAT: str = "{molecule_name}_GMX.gro"
    TOP_FILE_NAME_FORMAT: str = "{molecule_name}_GMX.top"
    POSRE_FILE_NAME_FORMAT: str = "posre_{molecule_name}.itp"
