import os
from enum import Enum
import dataclasses
from typing import Optional

<<<<<<< HEAD
TEMP_DIR = "temp"
=======
job_id = os.getenv("SLURM_JOB_ID", "local_run")
TEMP_DIR = f"temp_{job_id}"
>>>>>>> 91758eb (cleaned up)
LOG_DIR = "logs"
TOPOL_NAME = "topol.top"
MAIN_CACHE_DIR = "cache"
PREPROCESSED_DIR = "preprocessed"
PARAMETERISED_POLYMER_DIR = os.path.join(PREPROCESSED_DIR, "parameterised_polymers")
EQUILIBRIATED_SOLVENT_BOX_DIR = os.path.join(
    PREPROCESSED_DIR, "equilibriated_solvent_boxes"
)
SOLVENT_PDB_DIR = os.path.join(PREPROCESSED_DIR, "solvent_pdbs")
MDP_CACHE_DIR = os.path.join(MAIN_CACHE_DIR, "mdp_cache")
SHORT_POLYMER_BUILDING_BLOCKS_DIR = os.path.join(
    PREPROCESSED_DIR, "parameterised_polymer_building_blocks"
)

# SOLVENT_ITP_DIR = "preprocessed_output/solvents/solvent_itp"
# SOLVENT_MOL2_DIR = "preprocessed_output/solvents/solvent_mol2"
# SOLVENT_JSON_PATH = "preprocessed_output/solvents/solvent_atomtypes.json"


# ACPYPE_SOLVENT_OUTPUT_SUBDIR = "acpype_solvent_output"

ACPYPE_POLYMER_NAME = "POLY"
PACKMOL_TEMPLATE_DIR = "modules/packmol/templates"

EQUILIBRIATED_OUTPUTS_SUBDIR = "equilibriated_outputs"
MDP_TEMPLATE_DIR = "modules/gromacs/mdp_templates"
PME_TEMPLATE_DIR = os.path.join(MDP_TEMPLATE_DIR, "PME_mdps")
RF_TEMPLATE_DIR = os.path.join(MDP_TEMPLATE_DIR, "RF_mdps")
PREPROCESSED_PACKMOL_DIR = os.path.join(PREPROCESSED_DIR, "packmol")
