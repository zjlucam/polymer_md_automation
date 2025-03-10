from modules.workflows.base_workflow import BaseWorkflow
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from modules.cg_mappers.open_mscg_map_generator import OpenMSCGMapGenerator
from modules.utils.shared.file_utils import check_directory_exists, delete_directory
from modules.open_mscg.topol_generator import OpenMSCGTopolGenerator
from modules.open_mscg.trajectory_mapper import OpenMSCGTrajectoryMapper
from modules.open_mscg.force_matcher import OpenMSCGForceMatcher
from modules.open_mscg.topol_exporter import OpenMSCGTopolExporter
from config.paths import TEMP_DIR
from config.data_models.output_types import GromacsOutputs
from typing import List
import os


class CourseGrainer(BaseWorkflow):
    map_name = "cg_mapping"
    topol_name = "cg_topol"
    traj_name = "cg_traj"
    results_name = "results"
    cg_data_name = "cg_poly"
    include_ions_in_fm = False
    cg_non_bonding_padding_rule = "L2"
    cg_bonded_padding_rule = "LH"

    def __init__(
        self,
        polymer: BasePolymerGenerator,
        outputs: GromacsOutputs,
        output_dir: str,
        temp_dir: str = TEMP_DIR,
        polymer_resname: str = "UNL",
        ion_list: List[str] = ["NA", "CL"],
        sol_resname: str = "SOL",
        gro_to_open_mscg_offset: int = -1,
        model="BSpline",
        cleanup: bool = True,
        confirm_temp_dir_deletion: bool = True,
    ):
        super().__init__()
        self.polymer = polymer
        self.gro_file = outputs.gro
        self.trr_file = outputs.trr
        self.polymer_resname = polymer_resname
        self.ion_list = ion_list
        self.output_dir = output_dir
        self.model = model
        self.map_file = None
        self.cg_traj_file = None
        self.cg_topol_file = None
        self.cg_tables = None
        self.results = None
        self.mapper = None
        self.topol_generator = None
        self.traj_mapper = None
        self.force_matcher = None
        self.topol_exporter = None

        check_directory_exists(self.output_dir, make_dirs=True)
        self.subdir = self._generate_subdir()
        check_directory_exists(self.subdir, make_dirs=True)
        self.temp_dir = temp_dir
        check_directory_exists(self.temp_dir, make_dirs=True)
        self.sol_resname = sol_resname
        self.gro_to_open_mscg_offset = gro_to_open_mscg_offset

    def _generate_subdir(self):
        basename = os.path.basename(self.gro_file)
        subdir_name = basename.replace(".gro", "")
        self.subdir = os.path.join(self.output_dir, subdir_name)
        return self.subdir

    def run(self) -> str:
        self.traj_mapper = OpenMSCGTrajectoryMapper(
            polymer=self.polymer,
            gro_file=self.gro_file,
            polymer_resname=self.polymer_resname,
            ion_list=self.ion_list,
            sol_resname=self.sol_resname,
            gro_to_open_mscg_offset=self.gro_to_open_mscg_offset,
        )

        self.cg_traj_file = self.traj_mapper.run_cgmap(
            trr_path=self.trr_file,
            filename=self.traj_name,
            output_dir=self.temp_dir,
            map_filename=self.map_name,
        )
        self.map_file = self.traj_mapper.map_path
        self.topol_generator = OpenMSCGTopolGenerator(map_path=self.map_file)
        self.cg_topol_file = self.topol_generator.create_topol(
            filename=self.topol_name, output_dir=self.temp_dir
        )

        self.force_matcher = OpenMSCGForceMatcher(
            topol_generator=self.topol_generator,
            traj_path=self.cg_traj_file,
            sol_name=self.sol_resname,
            ions_list=self.ion_list,
            model=self.model,
        )

        results = self.force_matcher.run_cgfm(
            filename=self.results_name,
            output_dir=self.temp_dir,
            include_ions=self.include_ions_in_fm,
        )
        self.results = results
        tables_list = self.force_matcher.run_cgdump(
            output_dir=self.subdir,
            overwrite=True,
            non_bonding_padding_rule=self.cg_non_bonding_padding_rule,
            bonded_padding_rule=self.cg_bonded_padding_rule,
        )

        self.cg_tables = tables_list
        self.topol_exporter = OpenMSCGTopolExporter(
            cg_topol_file=self.cg_topol_file,
            cg_traj_file=self.cg_traj_file,
            cg_map_file=self.map_file,
        )

        topol_file = self.topol_exporter.run(
            filename=self.cg_data_name, output_dir=self.subdir
        )
        return topol_file

    def cleanup(self) -> None:
        if self.cleanup:
            delete_directory(
                TEMP_DIR, verbose=self.verbose, confirm=self.confirm_temp_dir_deletion
            )
