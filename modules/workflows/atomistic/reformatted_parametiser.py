from modules.rdkit.polymer_builders.alternating_copolymer import (
    AlternatingPolymerGenerator,
)
from config.acpype_config import AcpypeOutputConfig
from modules.file_conversion.converters.obabel_pdb_to_mol2_converter import (
    OBabelPDBtoMOL2Converter,
)
from modules.utils.shared.file_utils import (
    copy_file,
    delete_directory,
    check_directory_exists,
)
from typing import List, Optional
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.utils.atomistic.file_utils import calculate_minimum_box_size_from_df
from modules.gromacs.parsers.handlers.gro_handler import GroHandler
import re


from modules.acpype.acpype_parametizer import (
    ACPYPEParameterizer,
)
from config.data_models.output_types import GromacsPaths
from config.paths import (
    PARAMETERISED_POLYMER_DIR,
    TEMP_DIR,
    MAIN_CACHE_DIR,
    SHORT_POLYMER_BUILDING_BLOCKS_DIR,
)
from typing import List
import time as time
from abc import ABC, abstractmethod
from modules.workflows.base_workflow import BaseWorkflow
from modules.cache_store.pickle_cache import PickleCache
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from modules.rdkit.polymer_builders.alternating_copolymer import (
    AlternatingPolymerGenerator,
)
from modules.file_conversion.converters.editconf_pdb_to_gro import (
    EditconfPDBtoGROConverter,
)
from modules.rdkit.polymer_itp_scaler import PolymerITPScaler
import re
from typing import Dict, List, Union, Optional
from modules.gromacs.topol_generator.topol_generator import TopolGenerator
import logging
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# short_polymer_cache_dir = SHORT_POLYMER_CACHE_DIR
short_polymer_cache = PickleCache(name="short_polymer_cache", cache_dir=MAIN_CACHE_DIR)
# polymer_cache_dir = PARAMETERISED_POLYMER_DIR
long_polymer_cache_reformatted = PickleCache(
    cache_dir=MAIN_CACHE_DIR, name="long_polymer_cache_reformatted"
)


class ReformattedPolymerWorkflow(BaseWorkflow):
    box_dim_padding: float = 0.1

    def __init__(
        self,
        monomer_smiles: List[str],
        num_units: int,
        short_polymer_cache: PickleCache = short_polymer_cache,
        long_polymer_cache: PickleCache = long_polymer_cache_reformatted,
        polymer_generator: BasePolymerGenerator = AlternatingPolymerGenerator,
        verbose: bool = True,
        output_dir: str = PARAMETERISED_POLYMER_DIR,
        short_polymer_cache_dir: str = SHORT_POLYMER_BUILDING_BLOCKS_DIR,
        res_name="POLY",
    ):
        super().__init__()
        self.output_dir: str = output_dir
        self.short_polymer_cache_dir: str = short_polymer_cache_dir
        self.monomer_smiles: List[str] = monomer_smiles
        self.verbose: bool = verbose
        self.num_units: int = num_units
        self.short_polymer_cache = short_polymer_cache
        self.long_polymer_cache = long_polymer_cache
        self.res_name = res_name
        self.short_polymer_generator: BasePolymerGenerator = polymer_generator(
            monomer_smiles=self.monomer_smiles, res_name=self.res_name
        )
        self.long_polymer_generator: BasePolymerGenerator = polymer_generator(
            monomer_smiles=self.monomer_smiles, res_name=self.res_name
        )
        self.parametizer = ACPYPEParameterizer(acpype_molecule_name=self.res_name)
        self.file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
        self.num_repeats, self.actual_num_units = self._get_n_repeat()

    def check_short_polymer_cache(self):
        cache_key = "_".join(self.monomer_smiles)
        if self.short_polymer_cache.has_key(cache_key):
            parameterised_short_polymer = self.short_polymer_cache.retrieve_object(
                cache_key
            )
            logging.info(
                f"Short parameterised polymer retrieved from cache with key: {cache_key}"
            )
            return parameterised_short_polymer
        logging.info(
            f"Short parameterised polymer not found in cache with key: {cache_key}"
        )
        return None

    def _get_minimum_polymer_length(self):
        return len(self.monomer_smiles) + 2

    def _build_short_polymer(self, length: int):
        short_polymer_pdb = self.short_polymer_generator.generate_polymer(
            num_units=length,
            output_dir=TEMP_DIR,
            overwrite=True,
            save=True,
        )
        return short_polymer_pdb, self.short_polymer_generator.cg_map

    def build_and_parameterize_short_polymer(
        self,
        length: int,
        output_dir: str = SHORT_POLYMER_BUILDING_BLOCKS_DIR,
        mol2_output_dir: str = TEMP_DIR,
    ):
        pdb, short_cg_map = self._build_short_polymer(length)
        return self.parameterize_pdb(pdb, output_dir, mol2_output_dir), short_cg_map

    def parameterize_pdb(
        self, pdb_path: str, output_dir: str = TEMP_DIR, mol2_output_dir: str = TEMP_DIR
    ):
        check_directory_exists(output_dir)
        check_directory_exists(mol2_output_dir)
        converter = OBabelPDBtoMOL2Converter()
        mol2_file = converter.run(pdb_path, mol2_output_dir, verbose=self.verbose)
        return self.parametizer.run(
            mol2_file, output_dir, self.file_config, verbose=self.verbose
        )

    def _retrieve_or_build_short_parameterized_short_polymer(self):
        parameterised_short_polymer = self.check_short_polymer_cache()
        if parameterised_short_polymer:
            return (
                parameterised_short_polymer["outputs"],
                parameterised_short_polymer["cg_map"],
            )
        logging.info(f"Short polymer not found in cache, generating...")
        length = self._get_minimum_polymer_length()
        parameterised_files, short_cg_map = self.build_and_parameterize_short_polymer(
            length
        )

        cache_key = "_".join(self.monomer_smiles)
        self.short_polymer_cache.store_object(
            cache_key, {"outputs": parameterised_files, "cg_map": short_cg_map}
        )
        logger.info(f"Parameterised polymer saved to cache with key: {cache_key}")
        return parameterised_files, short_cg_map

    def _get_n_repeat(self):
        closest_multiple = (self.num_units - self._get_minimum_polymer_length()) // len(
            self.monomer_smiles
        )
        if closest_multiple < 1:
            length = self.num_units
        else:
            length = (
                len(self.monomer_smiles) * closest_multiple
                + self._get_minimum_polymer_length()
            )

        return closest_multiple, length

    def _extend_itp(
        self,
        output_path: str,
        short_polymer_itp: str,
        short_polymer_cg_map: Dict[str, List[Union[str, int]]],
        atom_start_index: Optional[int],
    ):

        itp_scaler = PolymerITPScaler(
            itp_path=short_polymer_itp,
            short_cg_map=short_polymer_cg_map,
            n_repeat=self.num_repeats,
            atom_start_index=atom_start_index,
        )
        itp_scaler.write_itp(output_path)
        logger.info(f"Long polymer built with {self.num_repeats} repeats")
        return output_path

    def _determine_atom_start_index(self, gro_path: str) -> Optional[int]:
        with open(gro_path, "r") as file:
            lines = file.readlines()

        first_atom_line = lines[2].strip().split()

        if len(first_atom_line) < 2:
            raise ValueError("Invalid .gro file format, atom name not found.")

        atom_name = first_atom_line[1]

        match = re.search(r"(\d+)$", atom_name)

        return int(match.group(1)) if match else None

    def _build_long_polymer(self, output_dir: str, pdb_output_dir: str = TEMP_DIR):
        short_polymer_files, cg_map = (
            self._retrieve_or_build_short_parameterized_short_polymer()
        )
        pdb = self.long_polymer_generator.generate_polymer(
            num_units=self.actual_num_units, output_dir=pdb_output_dir
        )
        pdb_basename = os.path.splitext(os.path.basename(pdb))[0]
        output_dir = os.path.join(output_dir, pdb_basename)
        gro = EditconfPDBtoGROConverter().run(pdb, output_dir)
        gro_path = os.path.join(output_dir, "POLY_GMX.gro")
        os.rename(gro, gro_path)
        gro = gro_path
        self.add_box_dim(gro_file=gro, padding=self.box_dim_padding)
        atom_start_index = self._determine_atom_start_index(gro)
        itp_output_path = self._get_itp_path(gro)
        itp = self._extend_itp(
            output_path=itp_output_path,
            short_polymer_itp=short_polymer_files.itp_path,
            short_polymer_cg_map=cg_map,
            atom_start_index=atom_start_index,
        )
        top_output_path = self._get_top_path(gro)
        top = TopolGenerator().create_topol(
            topol_path=top_output_path,
            res_name=self.res_name,
            itp_path=os.path.abspath(itp),
        )
        return GromacsPaths(
            itp_path=itp,
            gro_path=gro,
            top_path=top,
        )

    def _generate_polymer_cache_key(self, monomer_smiles: List[str], num_units: int):
        return "_".join(monomer_smiles) + "_" + str(num_units)

    def _get_itp_path(self, gro_path: str):
        return gro_path.replace(".gro", ".itp")

    def _get_top_path(self, gro_path: str):
        return gro_path.replace(".gro", ".top")

    def check_long_polymer_cache(self):
        cache_key = self._generate_polymer_cache_key(
            self.monomer_smiles, self.actual_num_units
        )
        if self.long_polymer_cache.has_key(cache_key):
            polymer, parameterised_polymer_pdb = (
                self.long_polymer_cache.retrieve_object(cache_key)
            )
            logging.info(
                f"Parameterised polymer retrieved from cache with key: {cache_key}"
            )
            return polymer, parameterised_polymer_pdb
        logging.info(f"Parameterised polymer not found in cache with key: {cache_key}")
        return None, None

    def run(self):
        output_dir = self.output_dir
        check_directory_exists(output_dir)
        polymer_generator, parameterised_polymer_pdb = self.check_long_polymer_cache()
        if polymer_generator and parameterised_polymer_pdb:
            self.long_polymer_generator = polymer_generator
            return parameterised_polymer_pdb
        logging.info(f"Long polymer not found in cache, generating...")

        if self.num_repeats < 1:
            polymer_paths, cg_map = self.build_and_parameterize_short_polymer(
                self.num_units, output_dir=output_dir
            )

        else:
            polymer_paths = self._build_long_polymer(output_dir)
        cache_key = self._generate_polymer_cache_key(
            self.monomer_smiles, self.actual_num_units
        )

        self.long_polymer_cache.store_object(
            cache_key, (self.long_polymer_generator, polymer_paths)
        )
        logger.info(f"Parameterised polymer saved to cache with key: {cache_key}")
        return polymer_generator

    @staticmethod
    def add_box_dim(
        gro_file, padding: float = 0.1, parser=GromacsParser(), gro_handler=GroHandler()
    ) -> List[float]:
        sections = parser.parse(gro_file)
        first_key = next(iter(sections))  # Get the first key
        gro_section = sections[first_key]
        gro_handler.process(gro_section)
        box_size = calculate_minimum_box_size_from_df(gro_handler.content, padding)
        gro_handler.box_dimensions = box_size
        sections[first_key] = gro_handler.export()

        parser.export(sections, gro_file)

        return gro_file
