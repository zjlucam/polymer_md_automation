from config.data_models.output_types import GromacsOutputs
from modules.gromacs.index_manager import GromacsIndexManager
from config.paths import TEMP_DIR
import subprocess
from typing import List
import os
import numpy as np
from scipy.signal import savgol_filter
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GromacsAnalyser:
    RADIUS_GYRATION_FILE = "gyrate.xvg"
    DIFFUSION_COEFFICIENT_FILE = "msd.xvg"
    SASA_FILE = "sasa.xvg"
    END_TO_END_DISTANCE_FILE = "end_to_end.xvg"

    def __init__(
        self,
        outputs: GromacsOutputs,
        poly_resname: str = "UNL",
        ion_resnames: List[str] = ["NA", "CL"],
        output_dir: str = TEMP_DIR,
    ):
        self.outputs = outputs
<<<<<<< HEAD
=======
        print(self.outputs)
>>>>>>> 91758eb (cleaned up)
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

        self.index_file = self._create_index_file(
            gro_file=outputs.gro,
            output_dir=self.output_dir,
            poly_resname=poly_resname,
            ion_resnames=ion_resnames,
        )

    def _create_index_file(
        self, gro_file: str, output_dir: str, poly_resname: str, ion_resnames: List[str]
    ) -> str:
        index_manager = GromacsIndexManager(
            gro_file=gro_file,
            output_dir=output_dir,
            poly_name=poly_resname,
            ion_names=ion_resnames,
        )
        return index_manager.create_index_file()

    def _get_output_path(self, basename: str) -> str:
        """
        Constructs the full output path inside output_dir.
        """
        return os.path.join(self.output_dir, basename)

    def extract_radius_of_gyration(self):
        output_path = self._get_output_path(self.RADIUS_GYRATION_FILE)

        cmd = [
            "gmx",
            "gyrate",
            "-f",
            self.outputs.xtc,
            "-s",
            self.outputs.tpr,
            "-n",
            self.index_file,
            "-o",
            output_path,
        ]

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate("Polymer\n")
        p.wait()

        data = np.loadtxt(output_path, comments=("#", "@"))
        Rg_mean = np.mean(data[:, 1])
        Rg_std = np.std(data[:, 1])
        return Rg_mean, Rg_std

    def extract_diffusion_coefficient(self):
        """
        Computes the diffusion coefficient D from mean squared displacement (MSD),
        applying a denoising step instead of error-raising.
        Uses:
        - Savitzky-Golay filter to smooth MSD fluctuations
        - Weighted least squares regression for a stable fit
        """
        output_path = self._get_output_path(self.DIFFUSION_COEFFICIENT_FILE)

        cmd = [
            "gmx",
            "msd",
            "-f",
            self.outputs.xtc,
            "-s",
            self.outputs.tpr,
            "-n",
            self.index_file,
            "-o",
            output_path,
        ]

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate("Polymer\n")
        p.wait()

        if not os.path.exists(output_path):
            logger.warning(f"MSD output file not found: {output_path}")
            return np.nan  # Return NaN if file is missing

        # Load MSD data
        try:
            data = np.loadtxt(output_path, comments=("#", "@"))
            t = data[:, 0]  # Time (ps)
            msd = data[:, 1]  # Mean Squared Displacement (nm²)
        except Exception as e:
            logger.warning(f"Error reading MSD file: {e}")
            return np.nan  # Return NaN if file is unreadable

        # **Step 1: Apply Savitzky-Golay filter for denoising**
        smoothed_msd = savgol_filter(msd, window_length=11, polyorder=2, mode="nearest")

        # **Step 2: Compute local slopes to find the stable diffusion regime**
        slopes = np.gradient(smoothed_msd, t)

        # Identify the region where slope is stable (avoid initial ballistic motion)
        stable_region = np.where((slopes > 0) & (slopes < 10 * np.mean(slopes)))[0]

        if len(stable_region) < 10:
            logger.warning(
                "Stable MSD region too short, using raw MSD data for fitting."
            )
            stable_region = np.arange(len(t))  # Use all data if filtering fails

        # **Step 3: Fit only within the stable diffusive regime**
        start_idx = stable_region[len(stable_region) // 3]
        end_idx = stable_region[-1]

        # Use weighted least squares to reduce the impact of noise
        weights = 1 / (1 + np.abs(t[start_idx:end_idx] - np.mean(t[start_idx:end_idx])))
        slope, intercept = np.polyfit(
            t[start_idx:end_idx], smoothed_msd[start_idx:end_idx], 1, w=weights
        )

        # Compute the diffusion coefficient (D = slope / 6 for 3D diffusion)
        D = slope / 6.0

        # **Step 4: Handle edge cases**
        if D < 0 or np.isnan(D):
            logger.warning(f"Invalid diffusion coefficient: {D} nm²/ps. Returning NaN.")
            return np.nan  # Instead of error, return NaN for robustness

        return D

    def extract_sasa(self):
        output_path = self._get_output_path(self.SASA_FILE)

        cmd = [
            "gmx",
            "sasa",
            "-f",
            self.outputs.xtc,
            "-s",
            self.outputs.tpr,
            "-n",
            self.index_file,
            "-o",
            output_path,
        ]

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate("Polymer\n")
        p.wait()

        data = np.loadtxt(output_path, comments=("#", "@"))
        SASA_mean = np.mean(data[:, 1])
        SASA_std = np.std(data[:, 1])
        return SASA_mean, SASA_std

    def extract_end_to_end_distance(self):
        output_path = self._get_output_path(self.END_TO_END_DISTANCE_FILE)

        cmd = [
            "gmx",
            "distance",
            "-f",
            self.outputs.xtc,
            "-s",
            self.outputs.tpr,
            "-n",
            self.index_file,
            "-oall",
            output_path,
        ]

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)

        selection_input = f"Polymer_Carbons_Start_and_End\n"
        p.communicate(input=selection_input)
        p.wait()

        if not os.path.exists(output_path):
<<<<<<< HEAD
            logger.warning(f"Output file not found: {output_path}")
            E2E_mean, E2E_std = 0, 0
=======
            raise FileNotFoundError(f"Output file not found: {output_path}")
>>>>>>> 91758eb (cleaned up)

        try:
            data = np.loadtxt(output_path, comments=("#", "@"))
            E2E_mean = np.mean(data[:, 1])
            E2E_std = np.std(data[:, 1])
            return E2E_mean, E2E_std
        except Exception as e:
<<<<<<< HEAD
            logger.error(f"Error reading end-to-end distance file: {e}")
            return 0, 0
=======
            raise RuntimeError(f"Error reading end-to-end distance file: {e}")
>>>>>>> 91758eb (cleaned up)
