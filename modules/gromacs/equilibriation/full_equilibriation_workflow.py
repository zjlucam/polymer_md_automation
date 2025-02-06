import shutil
import os
from typing import Dict, List, Optional
from modules.gromacs.equilibriation.base_workflow_step import BaseWorkflowStep
from modules.cache_store.mdp_cache import MDPCache
from modules.utils.shared.file_utils import (
    check_directory_exists,
    copy_and_rename,
)
from config.data_models.output_types import GromacsOutputs
from modules.utils.atomistic.mdp_utils import generate_dynamic_filename
from config.paths import EQUILIBRIATED_OUTPUTS_SUBDIR


class FullEquilibrationWorkflow:
    def __init__(self, mdp_cache: MDPCache):
        self.mdp_cache = mdp_cache
        self.em_steps = []  # Store EM steps separately
        self.thermal_steps = []  # Store temperature-dependent steps

    def add_em_step(
        self,
        step_name: str,
        workflow_step: BaseWorkflowStep,
        template_path: str,
        base_params: Dict[str, str],
<<<<<<< HEAD
        additional_flags = ["-nt", "6", "-ntomp", "6", "-pin", "on"]
    ):
        """Add energy minimization step."""
        self.em_steps.append((step_name, workflow_step, template_path, base_params, additional_flags))
=======
    ):
        """Add energy minimization step."""
        self.em_steps.append((step_name, workflow_step, template_path, base_params))
>>>>>>> 91758eb (cleaned up)

    def add_thermal_step(
        self,
        step_name: str,
        workflow_step: BaseWorkflowStep,
        template_path: str,
        base_params: Dict[str, str],
<<<<<<< HEAD
        additional_flags =  ["-nt", "8", "-ntomp", "8", "-pin", "on"]
    ):
        """Add thermal steps (e.g., NVT, NPT)."""
        self.thermal_steps.append(
            (step_name, workflow_step, template_path, base_params, additional_flags)
=======
    ):
        """Add thermal steps (e.g., NVT, NPT)."""
        self.thermal_steps.append(
            (step_name, workflow_step, template_path, base_params)
>>>>>>> 91758eb (cleaned up)
        )

    def run(
        self,
        input_gro_path: str,
        input_topol_path: str,
        temp_output_dir: str,
        main_output_dir: str,
        log_dir: str,
        varying_params_list: List[Dict[str, str]],
        files_to_keep: Optional[List[str]] = None,
        subdir: str = EQUILIBRIATED_OUTPUTS_SUBDIR,
        save_intermediate_edr: bool = True,
        save_intermediate_gro: bool = True,
        save_intermediate_log: bool = True,
        verbose: bool = True,
        file_name_override: Optional[str] = None,
    ):
        check_directory_exists(temp_output_dir)
        check_directory_exists(log_dir)
        check_directory_exists(main_output_dir)
        outputs = GromacsOutputs()
        if not files_to_keep:
            files_to_keep = ["gro"]

        main_output_dir = os.path.join(main_output_dir, subdir)
        os.makedirs(main_output_dir, exist_ok=True)

        current_gro_path = input_gro_path
        final_step_name = None

        # Run all EM steps first
<<<<<<< HEAD
        for step_name, step, template_path, base_params, additional_flags in self.em_steps:
=======
        for step_name, step, template_path, base_params in self.em_steps:
>>>>>>> 91758eb (cleaned up)
            current_gro_path = step.run(
                step_name=step_name,
                mdp_template_path=template_path,
                input_gro_path=current_gro_path,
                input_topol_path=input_topol_path,
                temp_output_dir=temp_output_dir,
                log_dir=log_dir,
                varying_params=base_params,  # No variation for EM steps
                mdp_cache=self.mdp_cache,
                save_intermediate_edr=save_intermediate_edr,
                save_intermediate_gro=save_intermediate_gro,
                save_intermediate_log=save_intermediate_log,
                verbose=verbose,
<<<<<<< HEAD
                additional_flags=additional_flags
=======
>>>>>>> 91758eb (cleaned up)
            )
            final_step_name = step_name  # Track the last step name

            if not current_gro_path:
                raise RuntimeError(
                    f"EM Step '{step_name}' did not generate a .gro file required for subsequent steps."
                )

        # Run thermal steps with varying parameters
        for varying_params in varying_params_list:
<<<<<<< HEAD
            for step_name, step, template_path, base_params, additional_flags in self.thermal_steps:
=======
            for step_name, step, template_path, base_params in self.thermal_steps:
>>>>>>> 91758eb (cleaned up)
                # Merge base and varying parameters
                params = {**base_params, **varying_params}

                # Run the step
                current_gro_path = step.run(
                    step_name=f"{step_name}",
                    mdp_template_path=template_path,
                    input_gro_path=current_gro_path,
                    input_topol_path=input_topol_path,
                    temp_output_dir=temp_output_dir,
                    log_dir=log_dir,
                    varying_params=params,
                    mdp_cache=self.mdp_cache,
                    save_intermediate_edr=save_intermediate_edr,
                    save_intermediate_gro=save_intermediate_gro,
                    save_intermediate_log=save_intermediate_log,
                    verbose=verbose,
<<<<<<< HEAD
                    additional_flags=additional_flags
=======
>>>>>>> 91758eb (cleaned up)
                )
                final_step_name = step_name  # Track the last step name

            if final_step_name and files_to_keep:
                for ext in files_to_keep:
                    file_name = f"{final_step_name}.{ext}"
                    file_path = os.path.join(temp_output_dir, file_name)
                    if file_name_override:
                        new_filename = file_name_override
                    else:
                        new_filename = generate_dynamic_filename(
                            varying_params, extension=None
                        )
                    if os.path.exists(file_path):
                        new_file_path = copy_and_rename(
                            file_path,
                            main_output_dir,
                            new_name=new_filename,
                            delete_original=False,
                            replace_if_exists=True,
                        )
                        if hasattr(outputs, ext):  # Ensure the extension is valid
                            setattr(outputs, ext, new_file_path)
                        else:
                            raise ValueError(
                                f"Extension '{ext}' is not a valid output type for GromacsOutputs."
                            )

        return main_output_dir, outputs
