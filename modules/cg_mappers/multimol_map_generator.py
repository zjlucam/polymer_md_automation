import yaml
import logging
import re
from typing import List, Dict, Optional, Tuple, Union
from rdkit import Chem
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.utils.shared.file_utils import check_directory_exists
from modules.cg_mappers.base_map_generator import BaseMapGenerator
from modules.utils.atomistic.file_utils import get_gro_handler
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from ruamel.yaml import YAML

logger = logging.getLogger(__name__)


class MultimolOpenMSCGMapGenerator(BaseMapGenerator):

    map_extension = "yaml"

    def __init__(
        self,
        polymer: BasePolymerGenerator,
        gro_file: str,
        polymer_resname: str = "UNL",
        ion_list: List[str] = ["NA", "CL"],
        sol_resname: str = "SOL",
        gro_to_open_mscg_offset: int = -1,
        sol_to_bead_ratio: int = 4,
    ):
        super().__init__(polymer)
        self.sol_to_bead_ratio = sol_to_bead_ratio
        self.gro_df = get_gro_handler(gro_file).content
        self.polymer_resname = polymer_resname
        self.ion_list = ion_list
        self.site_types = {}
        self.blocks: List[Dict[str, str]] = []
        self.sol_resname = sol_resname
        self.single_sol_resname = f"{sol_resname}_s"
        self._gro_to_open_mscg_offset = gro_to_open_mscg_offset

        self._process_gro()
        self._process_cg_map()

    def _validate_gro_headers(self):
        required_columns = {"Residue Name", "Residue Number", "Atom Name", "Atom Index"}
        if not required_columns.issubset(self.gro_df.columns):
            raise ValueError(
                f"Missing required columns in .gro file: {required_columns - set(self.gro_df.columns)}"
            )

    def _is_sol_molecule(self, resname: str) -> bool:
        if resname != self.polymer_resname and resname not in self.ion_list:
            return True
        else:
            return False

    def _add_site_type(self, bead_name: str, x_weight: List[float]) -> Dict:
        if bead_name not in self.site_types:
            self.site_types[bead_name] = self._generate_site_type(x_weight)

    def _generate_site_type(self, x_weight: List[float]) -> Dict:
        return {
            "index": list(range(len(x_weight))),
            "x-weight": [float(x_weight_val) for x_weight_val in x_weight],
            "f-weight": [1.0] * len(x_weight),
        }

    def _construct_sites(self, resname: str):
        if resname == self.polymer_resname:
            return [
                [cg_map["bead_type"], cg_map["atom_indices"][0]]
                for cg_map in self.bead_mappings
            ]
        else:
            return [[resname, 0]]

    @staticmethod
    def _is_new_molecule(prev_resnum: Optional[int], resnum: int) -> bool:
        return prev_resnum is not None and resnum != prev_resnum

    @staticmethod
    def _is_new_residue(prev_resname: str, resname: str) -> bool:
        return prev_resname != resname

    @staticmethod
    def _validate_x_weight(
        x_weight: List[float], previous_x_weight: Optional[List[float]]
    ) -> bool:
        if previous_x_weight is None:
            return True
        elif x_weight == previous_x_weight:
            return True
        else:
            raise ValueError("Inconsistent x-weight values in .gro file.")

    def _get_atomic_mass(self, atom_name: str) -> float:
        element = self.get_element(atom_name)
        return round(Chem.GetPeriodicTable().GetAtomicWeight(element))

    @staticmethod
    def get_element(atom_name: str):
        match = re.match(r"([A-Za-z]+)(\d*)", atom_name)
        element = match.group(1) if match else atom_name

        atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(element.capitalize())
        return Chem.GetPeriodicTable().GetElementSymbol(atomic_num)

    def _process_gro(self):
        (
            last_start_idx,
            last_resname,
            last_resnum,
            molecule_count,
            x_weight,
            previous_x_weight,
        ) = self._reset_residue_state(start_idx=0, resname=None)
        last_resname = None

        for resname, res_num, index, atom_name in zip(
            self.gro_df["Residue Name"],
            self.gro_df["Residue Number"],
            self.gro_df["Atom Index"],
            self.gro_df["Atom Name"],
        ):
            if self._is_new_molecule(last_resnum, res_num):
                if resname != self.polymer_resname:
                    self._validate_x_weight(
                        x_weight=x_weight, previous_x_weight=previous_x_weight
                    )

                molecule_count, x_weight, previous_x_weight = (
                    self._reset_molecule_state(
                        molecule_count=molecule_count, x_weight=x_weight
                    )
                )

            if self._is_new_residue(last_resname, resname):
                self._finalise_block(
                    last_resname=last_resname,
                    previous_x_weight=previous_x_weight,
                    start_index=last_start_idx,
                    molecule_count=molecule_count,
                )
                (
                    last_start_idx,
                    last_resname,
                    last_resnum,
                    molecule_count,
                    x_weight,
                    previous_x_weight,
                ) = self._reset_residue_state(start_idx=index, resname=resname)
            x_weight.append(self._get_atomic_mass(atom_name))
            last_resnum = res_num

        self._finalise_block(
            last_resname=last_resname,
            previous_x_weight=x_weight,
            start_index=last_start_idx,
            molecule_count=molecule_count + 1,
        )

    @staticmethod
    def _reset_residue_state(
        start_idx: int,
        resname: Optional[str] = None,
    ) -> Tuple[int, Optional[int], int, List[float], Optional[List[float]]]:
        last_resnum = None
        molecule_count = 0
        x_weight = []
        previous_x_weight = None
        last_resname = resname
        last_start_idx = start_idx
        return (
            last_start_idx,
            last_resname,
            last_resnum,
            molecule_count,
            x_weight,
            previous_x_weight,
        )

    @staticmethod
    def _reset_molecule_state(
        molecule_count: int,
        x_weight: List[float],
    ) -> Tuple[Optional[int], int, List[float], List[float]]:
        molecule_count += 1
        previous_x_weight = x_weight
        resetted_x_weight = []
        return molecule_count, resetted_x_weight, previous_x_weight

    def _finalise_block(
        self,
        last_resname: Optional[str],
        previous_x_weight: List[float],
        start_index: int,
        molecule_count: int,
    ) -> Dict[str, Union[str, List[float], int]]:
        if last_resname is None:
            return
        if self._is_sol_molecule(last_resname):
            resname_entry = self.single_sol_resname
        else:
            resname_entry = last_resname

        self.blocks.append(
            {
                "resname": resname_entry,
                "x-weight": previous_x_weight,
                "startidx": start_index
                + self._gro_to_open_mscg_offset,  # -1 as gro is 1-indexed, but openmscg is 0-indexed
                "molcount": molecule_count,
            }
        )
        if last_resname != self.polymer_resname:
            self._add_site_type(bead_name=resname_entry, x_weight=previous_x_weight)

    def _split_block(
        self, block: Dict[str, Union[str, List[float], int]]
    ) -> List[Dict[str, Union[str, List[float], int]]]:
        n_cg_mols, n_atom_mols = self._divmod(block["molcount"], self.sol_to_bead_ratio)
        cg_x_weight = block["x-weight"] * self.sol_to_bead_ratio
        cg_startidx = block["startidx"]

        cg_block = {
            "resname": self.sol_resname,
            "x-weight": cg_x_weight,
            "startidx": cg_startidx,
            "molcount": n_cg_mols,
        }
        self._add_site_type(bead_name=self.sol_resname, x_weight=cg_x_weight)

        if n_atom_mols == 0:
            return cg_block, None

        atom_startidx = len(cg_x_weight) + cg_startidx
        atom_block = {
            "resname": self.single_sol_resname,
            "x-weight": block["x-weight"],
            "startidx": atom_startidx,
            "molcount": n_atom_mols,
        }

        return cg_block, atom_block

    def _process_blocks_and_add_cg_sol(
        self, blocks: List[Dict[str, Union[str, List[float], int]]]
    ) -> List[Dict[str, Union[str, List[float], int]]]:
        logger.info(
            f"Adding CG solvent to mapping with {self.sol_to_bead_ratio} solvent molecules per bead."
        )
        i = 0
        while i < len(blocks):
            if blocks[i]["resname"] == self.single_sol_resname:
                cg_block, atom_block = self._split_block(blocks[i])
                blocks[i] = cg_block
                if atom_block is not None:
                    blocks.insert(i + 1, atom_block)
                    i += 1

            i += 1
        return blocks

    def _divmod(self, a, b):
        return a // b, a % b

    def _process_cg_map(self):
        for bead in self.bead_mappings:
            self._add_site_type(bead_name=bead["bead_type"], x_weight=bead["x-weight"])

    def _create_system_entry(self):
        if self.blocks is None:
            raise ValueError("No blocks to construct system entry.")

        if self.sol_to_bead_ratio > 1:
            self.blocks = self._process_blocks_and_add_cg_sol(self.blocks)

        system = []
        for block in self.blocks:
            system.append(
                {
                    "anchor": block["startidx"],
                    "repeat": block["molcount"],
                    "offset": len(block["x-weight"]),
                    "sites": self._construct_sites(block["resname"]),
                }
            )
        return system

    def _generate_mapping(self, start_index=None):
        mapping = {
            "site-types": self.site_types,
            "system": self._create_system_entry(),
        }

        yaml = YAML()
        yaml.default_flow_style = None  # Ensure structured formatting
        yaml.indent(mapping=2, sequence=4, offset=2)  # Control indentation

        # Dump YAML to a string
        from io import StringIO

        yaml_string = StringIO()
        yaml.dump(mapping, yaml_string)
        return yaml_string.getvalue()
