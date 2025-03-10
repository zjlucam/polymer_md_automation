from modules.cg_mappers.open_mscg_map_generator import OpenMSCGMapGenerator
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from modules.utils.shared.file_utils import check_directory_exists, check_file_type
from typing import List, Optional
from mscg import *
from mscg.cli import cgmap
import logging
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class OpenMSCGTrajectoryMapper:
    default_map_name = "open_mscg_map"
    default_traj_name = "open_mscg_traj"

    def __init__(
        self,
        polymer: BasePolymerGenerator,
        gro_file: str,
        polymer_resname: str = "UNL",
        ion_list: List[str] = ["NA", "CL"],
        sol_resname: str = "SOL",
        gro_to_open_mscg_offset: int = -1,
    ):
        self.map_path = None
        self.cg_lammpstrj = None
        self.map_generator: OpenMSCGMapGenerator = OpenMSCGMapGenerator(
            polymer=polymer,
            gro_file=gro_file,
            polymer_resname=polymer_resname,
            ion_list=ion_list,
            sol_resname=sol_resname,
            gro_to_open_mscg_offset=gro_to_open_mscg_offset,
        )

    def _generate_map(
        self, output_dir: Optional[str] = None, map_filename: Optional[str] = None
    ) -> str:
        if not map_filename:
            map_filename = self.default_map_name
        output_name = f"{map_filename}"
        map_path = self.map_generator.create_map(
            filename=output_name, output_dir=output_dir
        )
        return map_path

    def _generate_output_name(
        self,
        output_dir: Optional[str] = None,
        output_filename: Optional[str] = None,
    ):
        if not output_filename:
            output_filename = self.default_traj_name
        cg_traj_name = f"{output_filename}.lammpstrj"
        if not output_dir:
            traj_path = cg_traj_name
        else:
            traj_path = os.path.join(output_dir, cg_traj_name)

        return traj_path

    def run_cgmap(
        self,
        trr_path: str,
        filename: Optional[str] = None,
        output_dir: Optional[str] = None,
        map_filename: Optional[str] = None,
    ) -> str:
        check_file_type(trr_path, "trr")
        if output_dir:
            check_directory_exists(output_dir, make_dirs=True)
        map_path = self._generate_map(output_dir=output_dir, map_filename=map_filename)
        output_path = self._generate_output_name(
            output_dir=output_dir, output_filename=filename
        )
        self.map_path = map_path
        self.cg_lammpstrj = output_path

        logger.info(
            f"Running cgmap with map: {map_path}, traj: {trr_path}, out: {output_path}"
        )

        cgmap.main(map=map_path, traj=trr_path, out=output_path)
        return output_path

    def run_with_premade_map(
        self,
        trr_path: str,
        map_path: str,
        filename: Optional[str] = None,
        output_dir: Optional[str] = None,
    ) -> str:
        check_file_type(trr_path, "trr")
        if output_dir:
            check_directory_exists(output_dir, make_dirs=True)
        output_path = self._generate_output_name(
            output_dir=output_dir, output_filename=filename
        )
        self.map_path = map_path
        self.cg_lammpstrj = output_path

        logger.info(
            f"Running cgmap with map: {map_path}, traj: {trr_path}, out: {output_path}"
        )

        cgmap.main(map=map_path, traj=trr_path, out=output_path)
        return output_path
