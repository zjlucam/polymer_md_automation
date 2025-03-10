import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from modules.utils.shared.file_utils import prepare_output_file_path
from typing import Optional
import numpy as np
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def add_polymer_to_solvent(
    polymer_file: str,
    solvent_file: str,
    output_dir: Optional[str] = None,
    output_name: Optional[str] = None,
    cutoff=0.2,
) -> str:
    """
    Add a polymer to a pre-equilibrated solvent box, removing overlapping solvent molecules.

    :param polymer_file: Path to the polymer GRO file.
    :param solvent_file: Path to the solvent GRO file.
    :param output_dir: Directory to save the combined system.
    :param output_name: Name of the output file.
    :param cutoff: Distance cutoff (in nm) for removing overlapping solvent molecules.
    """
    output_path = prepare_output_file_path(polymer_file, "gro", output_dir, output_name)
    u_polymer = mda.Universe(polymer_file)
    u_solvent = mda.Universe(solvent_file)

    # Center polymer in the solvent box
    polymer_center = u_polymer.atoms.center_of_mass()
    solvent_center = u_solvent.dimensions[:3] / 2
    u_polymer.atoms.translate(solvent_center - polymer_center)

    # Identify overlapping solvent molecules
    distances = distance_array(
        u_polymer.atoms.positions, u_solvent.atoms.positions, box=u_solvent.dimensions
    )
    overlapping_atoms = distances.min(axis=0) < cutoff
    overlapping_residues = u_solvent.atoms[overlapping_atoms].residues.resids
    if overlapping_residues is None or len(overlapping_residues) == 0:
        logger.info("No overlapping residues detected, keeping all solvent molecules.")
        non_overlapping_solvent = u_solvent.atoms  # Keep all solvent molecules
    else:
        # Select non-overlapping solvent molecules
        non_overlapping_solvent = u_solvent.select_atoms(
            f"not resid {' '.join(map(str, overlapping_residues))}"
        )

    # Combine polymer and non-overlapping solvent
    combined = mda.Merge(u_polymer.atoms, non_overlapping_solvent)
    combined.dimensions = u_solvent.dimensions  # Retain original box dimensions

    # Save the combined system
    with mda.Writer(output_path, n_atoms=combined.atoms.n_atoms) as W:
        W.write(combined)

    return output_path


def add_n_parallel_polymers_to_solvent(
    polymer_file: str,
    solvent_file: str,
    n: int = 2,  # Number of polymers
    output_dir: Optional[str] = None,
    output_name: Optional[str] = None,
    cutoff=0.2,  # Solvent removal cutoff (nm)
    min_distance=1.5,  # Minimum distance between polymers (nm)
    align_axis="x",  # Primary alignment axis
) -> str:
    """
    Add `n` parallel polymer chains to a pre-equilibrated solvent box, ensuring no overlap.

    :param polymer_file: Path to the polymer GRO file.
    :param solvent_file: Path to the solvent GRO file.
    :param n: Number of polymer chains.
    :param output_dir: Directory to save the combined system.
    :param output_name: Name of the output file.
    :param cutoff: Distance cutoff (in nm) for removing overlapping solvent molecules.
    :param min_distance: Minimum separation distance between polymers.
    :param align_axis: Axis along which the polymers should be aligned ("x" or "z").
    :return: Path to the output `.gro` file.
    """
    output_path = prepare_output_file_path(polymer_file, "gro", output_dir, output_name)

    # Load solvent
    u_solvent = mda.Universe(solvent_file)
    box_length = u_solvent.dimensions[:3]
    solvent_center = box_length / 2

    # Load one polymer to get its dimensions
    u_polymer = mda.Universe(polymer_file)
    polymer_positions = u_polymer.atoms.positions
    polymer_width_min = np.min(polymer_positions, axis=0)
    polymer_width_max = np.max(polymer_positions, axis=0)
    polymer_size = polymer_width_max - polymer_width_min

    # Determine alignment axes
    primary_axis_idx = {"x": 0, "y": 1, "z": 2}[align_axis]
    offset_axis1 = (primary_axis_idx + 1) % 3  # First perpendicular axis
    offset_axis2 = (primary_axis_idx + 2) % 3  # Second perpendicular axis

    # Estimate grid layout (square-like)
    num_cols = int(np.ceil(np.sqrt(n)))  # Number of columns
    num_rows = int(np.ceil(n / num_cols))  # Number of rows

    # Compute polymer spacing
    polymer_spacing1 = polymer_size[offset_axis1] + min_distance
    polymer_spacing2 = polymer_size[offset_axis2] + min_distance

    # Create polymer instances
    polymers = []
    for i in range(n):
        polymer = mda.Universe(polymer_file)  # Create polymer instance

        # Compute row/column indices
        row_idx = i // num_cols
        col_idx = i % num_cols

        # Compute shift
        shift_vector = np.zeros(3)
        shift_vector[offset_axis1] = (row_idx - (num_rows / 2)) * polymer_spacing1
        shift_vector[offset_axis2] = (col_idx - (num_cols / 2)) * polymer_spacing2

        # Center polymer
        polymer_center = polymer.atoms.center_of_mass()
        polymer.atoms.translate(solvent_center - polymer_center + shift_vector)

        # Adjust residue IDs to prevent duplication
        max_resid = polymers[-1].residues.resids.max() if polymers else 0
        polymer.residues.resids += max_resid + 1

        polymers.append(polymer)

    # Ensure no polymer-polymer overlap
    for i in range(n):
        for j in range(i + 1, n):
            polymer_distance = np.min(
                distance_array(
                    polymers[i].atoms.positions,
                    polymers[j].atoms.positions,
                    box=u_solvent.dimensions,
                )
            )
            if polymer_distance < cutoff:
                logger.warning(f"⚠️ Polymers {i} and {j} are too close! Adjusting placement...")
                polymers[j].atoms.translate(min_distance)

    # Identify overlapping solvent molecules
    all_polymer_atoms = np.concatenate([p.atoms.positions for p in polymers])
    solvent_distances = distance_array(
        all_polymer_atoms, u_solvent.atoms.positions, box=u_solvent.dimensions
    )
    overlapping_atoms = solvent_distances.min(axis=0) < cutoff
    overlapping_residues = u_solvent.atoms[overlapping_atoms].residues.resids

    if overlapping_residues is None or len(overlapping_residues) == 0:
        logger.warning(" No overlapping residues detected, keeping all solvent molecules.")
        non_overlapping_solvent = u_solvent.atoms  # Keep all solvent molecules
    else:
        # Select non-overlapping solvent molecules
        non_overlapping_solvent = u_solvent.select_atoms(
            f"not resid {' '.join(map(str, overlapping_residues))}"
        )

    # Merge all polymers and the solvent
    combined = mda.Merge(*[p.atoms for p in polymers], non_overlapping_solvent)
    combined.dimensions = u_solvent.dimensions  # Retain original box dimensions

    # Save the combined system
    with mda.Writer(output_path, n_atoms=combined.atoms.n_atoms) as W:
        W.write(combined)

    return output_path
