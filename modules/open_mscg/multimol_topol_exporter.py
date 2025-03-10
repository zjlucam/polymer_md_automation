from mscg import Topology, Trajectory
from typing import List, Tuple, Optional
from modules.utils.shared.file_utils import check_directory_exists, move_files
import logging
import numpy as np
import yaml

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MultimolOpenMSCGTopolExporter:
    fmt = "lammpstrj"
    cg_system_name = "cg_poly"
    top_format = "lammps"

    def __init__(
        self,
        cg_topol_file: str,
        cg_traj_file: str,
        cg_map_file: str,
        overwrite: bool = True,
    ):
        self.map_bead_masses = None
        self.trajectory = None
        self.bead_mass_identities = None

        self.overwrite = overwrite

        self.cg_topol_file = cg_topol_file
        self._parse_top()
        self.cg_map_file = cg_map_file
        self._parse_map()
        self.cg_traj_file = cg_traj_file
        self._parse_traj()
        self._get_bead_masses()

    def _parse_map(self) -> None:
        with open(self.cg_map_file, "r") as file:
            data = yaml.safe_load(file)

        site_types = data.get("site-types", {})
        cg_bead_masses = [
            (site, sum(values.get("x-weight", [])))
            for site, values in site_types.items()
        ]
        self.map_bead_masses = cg_bead_masses

    def _parse_traj(self) -> Trajectory:
        CG_traj = Trajectory(self.cg_traj_file, fmt=self.fmt)
        CG_traj.read_frame()
        self.trajectory = CG_traj

    def _get_molecule_ids(
        self, molecule_count: int, num_beads_types: int
    ) -> List[List[int]]:
        upper_limit = molecule_count + 1
        return np.repeat(np.arange(1, upper_limit), num_beads_types).tolist()

    def _get_molecule_count(self, trajectory: Trajectory) -> int:
        return trajectory.x.shape[0]

    def _get_num_bead_types(self, bead_masses: List[Tuple[str, float]] = None) -> int:
        return len(bead_masses)

    def _create_topology(self, cg_system_name: Optional[str] = None) -> str:
        top = Topology.read_file(self.cg_topol_file)
        if not cg_system_name:
            cg_system_name = self.cg_system_name
        top.system_name = cg_system_name
        bead_masses = self._get_bead_masses()
        molecule_count = self._get_molecule_count(trajectory=self.trajectory)
        num_bead_types = self._get_num_bead_types(bead_masses=bead_masses)
        molecule_ids = self._get_molecule_ids(
            molecule_count=molecule_count, num_beads_types=num_bead_types
        )
        top.save(
            self.top_format,
            masses=bead_masses,
            molecule=molecule_ids,
            box=np.vstack([np.zeros(3), self.trajectory.box]).T,
            x=self.trajectory.x,
        )
        return f"{cg_system_name}.data"

    def _parse_top(self) -> None:
        with open(self.cg_topol_file, "r") as file:
            lines = file.readlines()

        cgtypes_section = False
        cg_order = []

        for line in lines:
            line = line.strip()
            if line.startswith("cgtypes"):
                cgtypes_section = True
                continue
            elif cgtypes_section and line.startswith("moltypes"):
                break
            elif cgtypes_section and line:
                cg_order.append(line)

        self.cg_order = cg_order

    def _get_bead_masses(self) -> List[float]:
        cg_order = self.cg_order

        cg_dict = {key: value for key, value in self.map_bead_masses}

        if set(cg_order) - set(cg_dict.keys()) != set():
            logger.warning("Mappings do not match between file and tuple.")
            logger.warning(f"Missing mappings: {set(cg_order) - set(cg_dict.keys())}")

        ordered_tuple = [(cg, cg_dict[cg]) for cg in cg_order if cg in cg_dict]
        self.bead_mass_identities = ordered_tuple
        masses = [mass for bead_name, mass in ordered_tuple]
        return masses

    def run(self, filename: str, output_dir: Optional[str] = None) -> str:
        top_path = self._create_topology(cg_system_name=filename)
        if output_dir:
            check_directory_exists(output_dir, make_dirs=True)
            top_path = move_files(
                file_paths=[top_path],
                target_directory=output_dir,
                overwrite=self.overwrite,
            )[0]
        return top_path
