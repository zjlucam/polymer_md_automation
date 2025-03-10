from typing import List, Optional
from modules.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from config.paths import TEMP_DIR, LOG_DIR, EQUILIBRIATED_OUTPUTS_SUBDIR
from modules.gromacs.commands.genion import GenIon
from modules.utils.shared.file_utils import copy_file
from config.data_models.output_types import GromacsPaths

import os
from config.data_models.output_types import GromacsOutputs
from config.mdp_workflow_config import minim_workflow, polymer_workflow
from modules.workflows.base_workflow import BaseWorkflow
from config.data_models.solvent import Solvent
from modules.workflows.atomistic.solvent_equilibriator import (
    SolventEquilibriationWorkflow,
)

from modules.utils.shared.file_utils import (
    prepare_output_file_path,
    check_directory_exists,
)
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from typing import Optional, Tuple, Dict
import pandas as pd
from modules.gromacs.parsers.handlers.data_handler import DataHandler
import logging
import os
import shutil
from config.data_models.output_types import GromacsPaths
from modules.utils.atomistic.file_utils import (
    get_gro_handler,
    calculate_molecule_counts,
    replace_dataframe_contents,
    replace_value_in_dataframe,
    create_includes_section,
    delete_all_include_sections,
    add_full_rows_to_handler_deduplicate,
)
from modules.workflows.atomistic.polymer_parametizer import PolymerGeneratorWorkflow
from modules.cache_store.equilibriated_atomistic_polymer_cache import (
    EquilibriatedAtomisticPolymerCache,
)
from modules.moltemplate.moltemplate_utils import (
    add_polymer_to_solvent,
)
from modules.utils.shared.file_utils import delete_directory
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
equiibriated_atomistic_polymer_cache = EquilibriatedAtomisticPolymerCache()


class PolymerEquilibriationWorkflow(BaseWorkflow):
    saved_file_types: List[str] = ["edr", "trr", "gro", "xtc", "tpr"]
    forcefield: str = "amber99sb-ildn.ff/forcefield.itp"
    ion_itp_file: str = "amber99sb-ildn.ff/ions.itp"
    topol_name: str = "topol"
    polymer_in_solvent_name: str = "polymer_in_solvent"
    polymer_addition_cutoff: float = 0.2

    def __init__(
        self,
        monomer_smiles: List[str],
        num_units: int,
        solvent: Solvent,
        box_size_nm: List[float],
        temperatures: List[float],
        output_dir: str,
        minim_workflow: FullEquilibrationWorkflow = minim_workflow,
        full_workflow: FullEquilibrationWorkflow = polymer_workflow,
        cache: EquilibriatedAtomisticPolymerCache = equiibriated_atomistic_polymer_cache,
        pos_ion_name: str = "NA",
        neg_ion_name: str = "CL",
        polymer_name: str = "POLY",
        verbose: bool = True,
        cleanup_log: bool = True,
        cleanup_temp: bool = True,
        confirm_temp_deletion: bool = True,
    ):
        super().__init__()
        self.verbose = verbose
        self.solvent = solvent
        self.outputs: List[GromacsOutputs] = []
        self.cache = cache
        self.polymer = None
        self.monomer_smiles = monomer_smiles
        self.num_units = num_units
        self.pname = pos_ion_name
        self.nname = neg_ion_name
        self.ions_list = [pos_ion_name, neg_ion_name]
        self.box_size_nm = box_size_nm
        self.temperatures = temperatures
        self.minim_workflow = minim_workflow
        self.full_workflow = full_workflow
        self.cleanup_log = cleanup_log
        self.cleanup_temp = cleanup_temp
        self.confirm_temp_deletion = confirm_temp_deletion
        self.polymer_name = polymer_name
        self.actual_num_units = None
        self.parameterised_polymer = self._retrieve_parameterised_polymer()
        self.subdir = self._retrieve_subdir_name(self.actual_num_units)
        self.final_output_dir = os.path.join(output_dir, self.subdir)
        check_directory_exists(TEMP_DIR, make_dirs=True)
        check_directory_exists(LOG_DIR, make_dirs=True)
        check_directory_exists(self.final_output_dir, make_dirs=True)

    def _retrieve_subdir_name(self, actual_num_units: int) -> str:
        monomer_smiles_str = "_".join(self.monomer_smiles)
        subdir = f"{self.solvent.name}_{self.solvent.compressibility}_{monomer_smiles_str}_{actual_num_units}"
        return subdir

    def _retrieve_parameterised_polymer(self) -> GromacsPaths:
        polymer_workflow = PolymerGeneratorWorkflow(
            monomer_smiles=self.monomer_smiles, num_units=self.num_units
        )

        self.actual_num_units = polymer_workflow.actual_num_units
        outputs = polymer_workflow.run()
        self.polymer = polymer_workflow.long_polymer_generator
        if outputs is None:
            raise ValueError("PolymerGeneratorWorkflow.run() returned None!")
        return outputs

    def check_polymer_cache(self, temperature: float) -> Optional[GromacsOutputs]:
        cache_key = self.cache.get_cache_key(
            solvent=self.solvent,
            monomer_smiles=self.monomer_smiles,
            num_units=self.num_units,
            temperature=temperature,
        )
        if self.cache.has_key(cache_key):
            parameterised_polymer = self.cache.retrieve_object(cache_key)
            logging.info(f"Polymer retrieved from cache with key: {cache_key}")
            return parameterised_polymer
        logging.info(f"Polymer not found in cache with key: {cache_key}")
        return None

    def _retrieve_solvent_box(self, temperature: float) -> GromacsOutputs:
        solvent_box = SolventEquilibriationWorkflow(
            solvent=self.solvent,
            box_size_nm=self.box_size_nm,
            temperatures=[temperature],
            confirm_temp_dir_deletion=self.confirm_temp_deletion,
        ).run()
        return solvent_box

    def run(self) -> List[GromacsOutputs]:
        for temperature in self.temperatures:
            outputs = self.check_polymer_cache(temperature)
            if outputs:
                self.outputs.append(outputs)
                continue
            logger.info(f"Polymer not found in cache, generating...")
            outputs = self._run_per_temp(temperature)
            cache_key = self.cache.get_cache_key(
                solvent=self.solvent,
                monomer_smiles=self.monomer_smiles,
                num_units=self.num_units,
                temperature=temperature,
            )
            self.outputs.append(outputs)
            self.cache.store_object(key=cache_key, data=outputs)

        if self.cleanup_log:
            delete_directory(LOG_DIR, verbose=self.verbose, confirm=False)
        if self.cleanup_temp:
            delete_directory(
                TEMP_DIR, verbose=self.verbose, confirm=self.confirm_temp_deletion
            )
        return self.outputs

    def _run_per_temp(self, temperature: float) -> GromacsOutputs:
        solvent_box = self._retrieve_solvent_box(temperature)
        parameterised_polymer = self._retrieve_parameterised_polymer()
        polymer_in_solvent = add_polymer_to_solvent(
            polymer_file=parameterised_polymer.gro_path,
            solvent_file=solvent_box.gro,
            output_dir=TEMP_DIR,
            output_name=self.polymer_in_solvent_name,
            cutoff=self.polymer_addition_cutoff,
        )
        initial_minim_files = self._prepare_solute_files(
            solute_itp_file=parameterised_polymer.itp_path,
            solvent_itp_file=solvent_box.itp,
            solvent_box_gro_file=polymer_in_solvent,
            input_top_file=parameterised_polymer.top_path,
            output_dir=TEMP_DIR,
        )

        print("!!!!!!!!5")
        print(initial_minim_files.gro_path)
        print("!!!!!!!!!!!!!5")
        _, outputs = self.minim_workflow.run(
            input_gro_path=initial_minim_files.gro_path,
            input_topol_path=initial_minim_files.top_path,
            main_output_dir=TEMP_DIR,
            temp_output_dir=TEMP_DIR,
            log_dir=LOG_DIR,
            varying_params_list=[None],
            files_to_keep=["gro", "tpr"],
            save_intermediate_edr=True,
            save_intermediate_gro=True,
            save_intermediate_log=True,
            subdir=EQUILIBRIATED_OUTPUTS_SUBDIR,
            verbose=self.verbose,
            file_name_override="initial_minim",
        )
        neutralised_gro = GenIon().run(
            input_box_gro_path=outputs.gro,
            tpr_path=outputs.tpr,
            top_path=initial_minim_files.top_path,
            pname=self.pname,
            nname=self.nname,
            output_dir=TEMP_DIR,
        )

        prepared_files = self._prepare_solute_files(
            solute_itp_file=self.parameterised_polymer.itp_path,
            solvent_itp_file=solvent_box.itp,
            solvent_box_gro_file=neutralised_gro,
            input_top_file=self.parameterised_polymer.top_path,
            output_dir=TEMP_DIR,
        )
        print("!!!!!!!6")
        print(prepared_files)
        varying_params_list = self._create_varying_params_list(temperature)
        _, outputs = self.full_workflow.run(
            input_gro_path=prepared_files.gro_path,
            input_topol_path=prepared_files.top_path,
            temp_output_dir=TEMP_DIR,
            main_output_dir=self.final_output_dir,
            log_dir=LOG_DIR,
            files_to_keep=self.saved_file_types,
            save_intermediate_edr=True,
            save_intermediate_gro=True,
            save_intermediate_log=True,
            subdir=EQUILIBRIATED_OUTPUTS_SUBDIR,
            verbose=self.verbose,
            varying_params_list=varying_params_list,
        )
        topol_file = copy_file(
            initial_minim_files.top_path, self.final_output_dir, skip_if_exists=True
        )
        outputs.top = topol_file

        return outputs

    def _create_varying_params_list(self, temperature: float) -> List[Dict[str, str]]:
        return [
            {
                "temp": str(temperature),
                "compressibility": str(self.solvent.compressibility),
            }
        ]

    def _prepare_solute_files(
        self,
        solute_itp_file: str,
        solvent_itp_file: str,
        solvent_box_gro_file: str,
        input_top_file: str,
        output_dir: str,
        output_solute_itp_name: Optional[str] = None,
        output_solvent_itp_name: Optional[str] = None,
        parser: GromacsParser = GromacsParser(),
    ) -> GromacsPaths:
        parser = GromacsParser()
        solute_molecule_name = self.polymer_name

        gro_handler = get_gro_handler(solvent_box_gro_file)

        molecule_content = calculate_molecule_counts(
            gro_handler=gro_handler,
            residue_name_col="Residue Name",
            residue_number_col="Residue Number",
        )
        molecule_content = replace_value_in_dataframe(
            molecule_content,
            target_value="UNL",
            replacement_value=solute_molecule_name,
            move_to_top=True,
        )

        solute_itp_file, solvent_itp_file = self._process_solute_and_solvent_itps(
            solute_itp_file=solute_itp_file,
            solvent_itp_file=solvent_itp_file,
            output_dir=output_dir,
            solute_output_name=output_solute_itp_name,
            solvent_output_name=output_solvent_itp_name,
            parser=parser,
        )

        solute_itp_path = os.path.abspath(solute_itp_file)
        solvent_itp_path = os.path.abspath(solvent_itp_file)

        output_top_path = self._prepare_solute_topol(
            input_top_file=input_top_file,
            new_molecule_dataframe=molecule_content,
            forcefield=self.forcefield,
            ions_itp_file=self.ion_itp_file,
            solute_itp_file=solute_itp_path,
            solvent_itp_file=solvent_itp_path,
            output_name=self.topol_name,
            output_dir=output_dir,
            parser=parser,
            del_posre=True,
            del_defaults=True,
        )

        output_gro_path = prepare_output_file_path(
            solvent_box_gro_file, output_extension="gro", output_dir=output_dir
        )
        if output_gro_path != solvent_box_gro_file:
            shutil.copy(solvent_box_gro_file, output_gro_path)
        paths = GromacsPaths(
            itp_path=solute_itp_file, gro_path=output_gro_path, top_path=output_top_path
        )
        return paths

    def _process_solvent_itp(
        self,
        solvent_itp_file: str,
        output_dir: Optional[str] = None,
        output_name: Optional[str] = None,
        parser: GromacsParser = GromacsParser(),
        data_handler: DataHandler = DataHandler,
    ) -> Tuple[str, pd.DataFrame]:
        data_handler = data_handler()

        solvent_sections = parser.parse(solvent_itp_file)
        atoms_sections = solvent_sections["data_atomtypes"]
        data_handler.process(atoms_sections)
        atom_content = data_handler.content
        solvent_sections.pop("data_atomtypes")
        output_itp_path = prepare_output_file_path(
            solvent_itp_file, "itp", output_dir, output_name
        )
        output_itp_path = parser.export(solvent_sections, output_itp_path)
        return output_itp_path, atom_content

    def _process_solute_itp(
        self,
        solute_itp_file: str,
        solvent_atomtype_data: pd.DataFrame,
        output_dir: Optional[str] = None,
        output_name: Optional[str] = None,
        parser: GromacsParser = GromacsParser(),
        data_handler: DataHandler = DataHandler,
    ) -> str:
        data_handler = data_handler()
        solute_sections = parser.parse(solute_itp_file)
        atoms_sections = solute_sections["data_atomtypes"]

        data_handler.process(atoms_sections)
        data_handler = add_full_rows_to_handler_deduplicate(
            data_handler,
            solvent_atomtype_data,
            add_to_top=False,
            deduplicate_column="name",
        )
        atoms_sections = data_handler.export()
        solute_sections["data_atomtypes"] = atoms_sections

        output_itp_path = prepare_output_file_path(
            solute_itp_file, "itp", output_dir, output_name
        )
        output_itp_path = parser.export(solute_sections, output_itp_path)
        return output_itp_path

    def _process_solute_and_solvent_itps(
        self,
        solute_itp_file: str,
        solvent_itp_file: str,
        output_dir: str,
        solute_output_name: str,
        solvent_output_name: str,
        parser: GromacsParser = GromacsParser(),
        data_handler: DataHandler = DataHandler,
    ) -> Tuple[str, str]:
        output_solvent_itp, solvent_atom_data = self._process_solvent_itp(
            solvent_itp_file, output_dir, solvent_output_name, parser=parser
        )
        output_solute_itp = self._process_solute_itp(
            solute_itp_file,
            solvent_atom_data,
            output_dir=output_dir,
            output_name=solute_output_name,
            parser=parser,
        )

        return output_solute_itp, output_solvent_itp

    def _prepare_solute_topol(
        self,
        input_top_file: str,
        new_molecule_dataframe: pd.DataFrame,
        forcefield: str,
        ions_itp_file: str,
        solute_itp_file: str,
        solvent_itp_file: str,
        output_name: Optional[str] = None,
        output_dir: Optional[str] = None,
        parser: GromacsParser = GromacsParser(),
        del_posre: bool = True,
        del_defaults: bool = True,
    ):

        sections = parser.parse(input_top_file)

        sections = delete_all_include_sections(sections)

        sections["include_ions_itp"] = create_includes_section(ions_itp_file)
        sections.move_to_end("include_ions_itp", last=False)

        sections["include_solvent_itp"] = create_includes_section(solvent_itp_file)
        sections.move_to_end("include_solvent_itp", last=False)

        sections["include_solute_itp"] = create_includes_section(solute_itp_file)
        sections.move_to_end("include_solute_itp", last=False)

        sections["include_forcefield"] = create_includes_section(forcefield)
        sections.move_to_end("include_forcefield", last=False)

        if "data_molecules" not in sections:
            raise ValueError("No 'data_molecules' section found in topology file.")
        else:
            data_molecules_section = sections["data_molecules"]
            data_molecules_handler = parser.handler_registry.get_handler(
                data_molecules_section.construct_name
            )()
            data_molecules_handler.process(data_molecules_section)
            new_content = replace_dataframe_contents(
                original_df=data_molecules_handler.content,
                new_df=new_molecule_dataframe,
                pad_missing=True,
            )
            data_molecules_handler.content = new_content
            data_molecules_section = data_molecules_handler.export()
            sections["data_molecules"] = data_molecules_section

        # Optionally delete `data_posre` section
        if del_posre and "conditional_if" in sections:
            del sections["conditional_if"]

        if del_defaults and "data_defaults" in sections:
            del sections["data_defaults"]

        output_path = prepare_output_file_path(
            input_top_file, "top", output_dir, output_name
        )
        output_path = parser.export(sections, output_path)

        return output_path
