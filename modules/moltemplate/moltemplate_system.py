import math
import numpy as np
import logging
from typing import List
from scipy.spatial import cKDTree
from modules.moltemplate.polymer import MoltemplatePolymer
from modules.moltemplate.solvent import MoltemplateSolvent
from modules.moltemplate.multimol_solvent import MultimolMoltemplateSolvent
from config.data_models.solvent import Solvent
from modules.utils.shared.file_utils import check_directory_exists, check_file_type
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from modules.lammps.parsers.open_mscg_data_parser import OpenMSCGDataParser
from config.constants import ANGSTROM3_TO_M3, AVOGADROS_NUMBER, KG_TO_G

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MoltemplateSystem:
    units = "angstrom"

    def __init__(
        self,
        n_units: int,
        polymer: BasePolymerGenerator,
        box_nm: List[float],
        solvent: Solvent,
        openmscg_topol_path: str,
        default_bond_length: float = 3.0,
        padding: float = 1.5,  # Overlap cutoff (Å)
        sol_resname: str = "SOL",
        min_box_factor: float = 1.2,
        auto_gen_valid_dims: bool = True,
        sol_to_bead_ratio: int = 4,
    ):
        check_file_type(openmscg_topol_path, "data")
        self.sol_to_bead_ratio = sol_to_bead_ratio
        self.padding = padding
        self.min_box_factor = min_box_factor
        self.solvent: Solvent = solvent
        self.openmscg_parser = OpenMSCGDataParser(
            data_file=openmscg_topol_path, default_bond_length=default_bond_length
        )

        self.polymer_obj = MoltemplatePolymer(
            n_units=n_units,
            polymer_generator=polymer,
            openmscg_parser=self.openmscg_parser,
        )

        self.solvent_obj = MultimolMoltemplateSolvent(
            solvent=solvent,
            openmscg_parser=self.openmscg_parser,
            sol_resname=sol_resname,
            sol_to_bead_ratio=sol_to_bead_ratio,
        )

        self.box_nm = box_nm
        self.box_angstroms = self._validate_and_convert_box_dims(
            auto_gen_valid_dims=auto_gen_valid_dims
        )

        self.solvent_positions = []
        self.polymer_position = self._calculate_box_center()
        self.polymer_coords = self._extract_polymer_coords()
        self.polymer_tree = cKDTree(self.polymer_coords)

    def _calculate_box_center(self) -> tuple:
        """
        Computes the center of the box in Angstroms.
        """
        width, length, height = self.box_angstroms
        return (width / 2, length / 2, height / 2)

    def _extract_polymer_coords(self) -> list:
        """
        Extracts polymer bead positions dynamically from `lt_block`,
        ensuring that we correctly parse positions rather than assuming
        fixed bond lengths.
        """
        lt_lines = self.polymer_obj.lt_block.split("\n")

        coords = []
        for line in lt_lines:
            parts = line.strip().split()
            if parts and parts[0].startswith("$atom:"):
                try:
                    x, y, z = map(float, parts[1:4])
                    coords.append((x, y, z))
                except ValueError:
                    continue  # Skip invalid lines

        if not coords:
            raise ValueError("No valid polymer coordinates were found in `lt_block`.")

        # Calculate center offset for correct centering
        polymer_min = np.min(coords, axis=0)
        polymer_max = np.max(coords, axis=0)
        polymer_center = (polymer_min + polymer_max) / 2

        # Shift polymer to center of the box
        box_center = np.array(self.polymer_position)
        centered_coords = [
            tuple(coord - polymer_center + box_center) for coord in coords
        ]

        return centered_coords

    def _get_average_bond_length(self):
        bond_lengths = self.openmscg_parser.bond_lengths
        avg_bond_length = np.mean(list(bond_lengths.values())[1:-1])
        return avg_bond_length

    def _validate_and_convert_box_dims(
        self, auto_gen_valid_dims: bool = True
    ) -> List[float]:
        box_angstroms = self.convert_box_dim_nm_to_ang(self.box_nm)
        min_box_dim = (
            self.polymer_obj.actual_n
            * self._get_average_bond_length()
            * self.min_box_factor
        )
        if any(dim < min_box_dim for dim in box_angstroms):
            logger.warning(
                f"Box dimensions {self.box_nm[0]:.1f}nm x {self.box_nm[1]:.1f}nm x {self.box_nm[2]:.1f}nm are too small. Minimum box dimension: {min_box_dim:.1f}Å or {(min_box_dim/10):.1f}nm"
            )
            if auto_gen_valid_dims:
                logger.info("Attempting to generate valid box dimensions...")
                box_angstroms = self._generate_valid_box_dims(
                    box_dims_angstrom=box_angstroms, min_box_dim=min_box_dim
                )
                logger.info(
                    f"New box dimensions:{box_angstroms[0]*0.1:.1f}nm x {box_angstroms[1]*0.1:.1f}nm x {box_angstroms[2]*0.1:.1f}nm "
                )
                return box_angstroms
            else:
                raise ValueError("Please provide a larger box size.")

    @staticmethod
    def _generate_valid_box_dims(
        box_dims_angstrom: List[float], min_box_dim: float
    ) -> List[float]:
        for i, dim in enumerate(box_dims_angstrom):
            if dim < min_box_dim:
                box_dims_angstrom[i] = min_box_dim
        return box_dims_angstrom

    @staticmethod
    def convert_box_dim_nm_to_ang(box_nm: List[float]) -> None:
        return [dim * 10 for dim in box_nm]

    def _place_solvent(self):
        n_sol = self._calculate_num_solvent_molecules()
        ###########################################################
        positions = self._generate_solvent_positions_grid(
            n_sol=n_sol, d_cut=50, beta=0.05, min_density=0.05
        )
        ###########################################################
        polymer_coords = self.polymer_coords
        filtered_positions = self._remove_overlaps(
            sol_positions=positions, poly_positions=polymer_coords, cutoff=self.padding
        )

        self.solvent_positions = filtered_positions
        logging.info(
            f"Final solvent count after overlap removal: {len(self.solvent_positions)}"
        )

    def _calculate_num_solvent_molecules(self) -> int:

        width_ang, length_ang, height_ang = self.box_angstroms
        volume_ang = width_ang * length_ang * height_ang
        volume_m3 = volume_ang * ANGSTROM3_TO_M3

        mass_kg = self.solvent.density * volume_m3
        mass_g = mass_kg * KG_TO_G
        moles = mass_g / self.solvent_obj.mass
        n_molecules = int(math.floor(moles * AVOGADROS_NUMBER))

        logging.info(
            f"Calculated solvent molecules required for box {self.box_angstroms[0]*10:.1f}nm x {self.box_angstroms[1]*10:.1f}nm x {self.box_angstroms[2]*10:.1f}nm: {n_molecules}"
        )
        return n_molecules

    def _generate_solvent_positions_grid(
        self, n_sol: int, d_cut: float, beta: float, min_density: float
    ) -> list:
        """
        Generates solvent positions using a hybrid density approach:
        - Full density within `d_cut` (default 3 nm).
        - Rapid density decay after `d_cut`.
        - Uses exponential decay beyond cutoff.

        Parameters:
            - d_cut: Distance (Å) within which full density is maintained.
            - beta: Decay rate (higher = faster drop-off).
            - min_density: Minimum density floor (prevents depletion).

        Returns:
            - List of solvent positions.
        """

        logging.info("Generating solvent positions using hybrid density approach...")

        w, l, h = self.box_angstroms

        # Default uniform spacing estimate
        spacing = (w * l * h / n_sol) ** (1 / 3)
        nx, ny, nz = int(w // spacing), int(l // spacing), int(h // spacing)

        x0, x1 = -0.5 * w, 0.5 * w
        y0, y1 = -0.5 * l, 0.5 * l
        z0, z1 = -0.5 * h, 0.5 * h

        # Generate an initial uniform solvent grid
        all_positions = [
            (xx, yy, zz)
            for xx in np.linspace(x0, x1, nx)
            for yy in np.linspace(y0, y1, ny)
            for zz in np.linspace(z0, z1, nz)
        ]

        # KD-Tree for efficient polymer distance checking
        polymer_tree = self.polymer_tree

        # Define hybrid solvent density function
        def solvent_density_factor(pos):
            distance = polymer_tree.query(pos)[0]  # Find nearest polymer distance
            if distance <= d_cut:
                return 1  # Full density inside cutoff
            else:
                return max(np.exp(-beta * (distance - d_cut)), min_density)

        # Filter solvent positions adaptively
        final_positions = [
            pos
            for i, pos in enumerate(all_positions)
            if np.random.rand() < solvent_density_factor(pos)
        ]

        logging.info(
            f"Generated {len(final_positions)} solvent positions after hybrid adaptive filtering."
        )

        return final_positions

    def _remove_overlaps(
        self, sol_positions: list, poly_positions: list, cutoff: float
    ) -> list:

        logger.info("Running optimized solvent overlap removal...")

        # Convert lists to numpy arrays for efficiency
        sol_positions_np = np.array(sol_positions)
        poly_positions_np = np.array(poly_positions)

        # Build a KD-Tree for polymer positions
        polymer_tree = cKDTree(poly_positions_np)

        # Query tree: Find all solvent molecules within cutoff distance of any polymer bead
        close_solvent_indices = polymer_tree.query_ball_point(sol_positions_np, cutoff)

        # Flatten list and get unique solvent indices to remove
        close_solvent_indices = set(
            idx for sublist in close_solvent_indices for idx in sublist
        )

        # Keep only solvents **not** in the overlap list
        filtered_positions = [
            sol_positions[i]
            for i in range(len(sol_positions))
            if i not in close_solvent_indices
        ]

        logger.info(
            f"Removed {len(sol_positions) - len(filtered_positions)} overlapping solvents."
        )
        return filtered_positions

    def write_system_lt(self, filename="system.lt", output_dir=None):
        if output_dir:
            check_directory_exists(output_dir, make_dirs=True)
            filename = f"{output_dir}/{filename}"

        self._place_solvent()
        with open(filename, "w") as f:
            f.write("# Moltemplate system file\n")
            f.write("# Polymer is centered in the box\n")
            f.write("# Interactions are defined via tabulated .table files\n\n")

            # place polymer in center of box
            f.write(self.polymer_obj.lt_block)
            f.write(
                f"new Polymer polymer_instance {{ translate {self.polymer_position[0]:.2f} {self.polymer_position[1]:.2f} {self.polymer_position[2]:.2f} }}\n\n"
            )

            # solvents placed in b ox
            if self.solvent:
                f.write(self.solvent_obj.lt_block)
                for i, (sx, sy, sz) in enumerate(self.solvent_positions):
                    f.write(
                        f"new Solvent solvent_{i+1} {{ translate {sx:.2f} {sy:.2f} {sz:.2f} }}\n"
                    )

        logging.info(f"Moltemplate system file '{filename}' generated.")
