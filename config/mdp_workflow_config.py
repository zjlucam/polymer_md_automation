from config.paths import MDP_CACHE_DIR, RF_TEMPLATE_DIR, PME_TEMPLATE_DIR
from modules.gromacs.equilibriation.base_workflow_step import (
    BaseWorkflowStep,
)
from modules.gromacs.commands.grompp import Grompp
from modules.gromacs.commands.mdrun import MDrun
from modules.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from modules.cache_store.mdp_cache import MDPCache
import os

mdp_cache = MDPCache(cache_dir=MDP_CACHE_DIR)
workflow_step = BaseWorkflowStep(Grompp(), MDrun())
solvent_workflow = FullEquilibrationWorkflow(mdp_cache)


solvent_workflow.add_em_step(
    step_name="soft_em",
    workflow_step=workflow_step,
    template_path=os.path.join(PME_TEMPLATE_DIR, "em.mdp"),
    base_params={"nsteps": "10000", "emtol": "1000"},
)
solvent_workflow.add_em_step(
    step_name="minim_1",
    workflow_step=workflow_step,
    template_path=os.path.join(PME_TEMPLATE_DIR, "em.mdp"),
    base_params={"nsteps": "50000", "emtol": "1000"},
)
solvent_workflow.add_em_step(
    step_name="minim_2",
    workflow_step=workflow_step,
    template_path=os.path.join(PME_TEMPLATE_DIR, "em.mdp"),
    base_params={"nsteps": "20000", "emtol": "100"},
)

solvent_workflow.add_thermal_step(
    step_name="nvt",
    workflow_step=workflow_step,
    template_path=os.path.join(PME_TEMPLATE_DIR, "nvt.mdp"),
    base_params={
        "nsteps": "50000",
    },
)

solvent_workflow.add_thermal_step(
    step_name="npt_c_rescale_short",
    workflow_step=workflow_step,
    template_path=os.path.join(PME_TEMPLATE_DIR, "npt_c_rescale.mdp"),
    base_params={
        "nsteps": "20000",
        "dt": "0.002",
    },
)
solvent_workflow.add_thermal_step(
    step_name="npt_PR_long",
    workflow_step=workflow_step,
    template_path=os.path.join(PME_TEMPLATE_DIR, "npt_PR.mdp"),
    base_params={
        "nsteps": "15000",
        "dt": "0.001",
    },
)

minim_workflow = FullEquilibrationWorkflow(mdp_cache)
minim_workflow.add_em_step(
    step_name="em_init",
    workflow_step=workflow_step,
    template_path=os.path.join(RF_TEMPLATE_DIR, "em.mdp"),
    base_params={"nsteps": "50000", "emtol": "1000"},
)


polymer_workflow = FullEquilibrationWorkflow(mdp_cache)


polymer_workflow.add_em_step(
    step_name="minim_1_RF",
    workflow_step=workflow_step,
    template_path=os.path.join(RF_TEMPLATE_DIR, "em.mdp"),
    base_params={"nsteps": "50000", "emtol": "1000"},
)
polymer_workflow.add_em_step(
    step_name="minim_2_RF",
    workflow_step=workflow_step,
    template_path=os.path.join(RF_TEMPLATE_DIR, "em.mdp"),
    base_params={"nsteps": "20000", "emtol": "100"},
)

polymer_workflow.add_thermal_step(
    step_name="nvt_RF",
    workflow_step=workflow_step,
    template_path=os.path.join(RF_TEMPLATE_DIR, "nvt.mdp"),
    base_params={
        "nsteps": "30000",
    },
)

polymer_workflow.add_thermal_step(
    step_name="npt_c_rescale_short_RF",
    workflow_step=workflow_step,
    template_path=os.path.join(RF_TEMPLATE_DIR, "npt_c_rescale.mdp"),
    base_params={
        "nsteps": "20000",
        "dt": "0.002",
    },
)
polymer_workflow.add_thermal_step(
    step_name="npt_PR_long_RF",
    workflow_step=workflow_step,
    template_path=os.path.join(RF_TEMPLATE_DIR, "npt_PR.mdp"),
    base_params={
        "nsteps": "15000",
        "dt": "0.001",
    },
)
polymer_workflow.add_thermal_step(
    step_name="production_RF",
    workflow_step=workflow_step,
    template_path=os.path.join(RF_TEMPLATE_DIR, "prod.mdp"),
    base_params={
        "nsteps": "80000",
        "dt": "0.002",
    },
)
