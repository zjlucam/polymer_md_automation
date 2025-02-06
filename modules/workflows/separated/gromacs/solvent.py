from modules.utils.shared.file_utils import (
    check_file_type,
    prepare_output_file_path,
    delete_directory,
    check_directory_exists,
)
from modules.file_conversion.converters.obabel_pdb_to_mol2_converter import (
    OBabelPDBtoMOL2Converter,
)
from config.acpype_config import AcpypeOutputConfig
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from config.data_models.solvent import Solvent
from config.paths import (
    EQUILIBRIATED_SOLVENT_BOX_DIR,
    EQUILIBRIATED_OUTPUTS_SUBDIR,
    TEMP_DIR,
    LOG_DIR,
    PREPROCESSED_PACKMOL_DIR,
)
from typing import Optional
import logging
import os
from config.data_models.output_types import GromacsPaths, GromacsOutputs
from modules.utils.atomistic.file_utils import (
    get_gro_handler,
    get_residue_number,
    rename_residue_name_from_handler,
    rename_data_column_content,
    export_gro_handler,
    create_includes_section,
    delete_all_include_sections,
)
from modules.cache_store.pickle_cache import PickleCache
from modules.file_conversion.converters.base_converter import BaseConverter
from modules.file_conversion.converters.editconf_gro_to_pdb import (
    EditconfGROtoPDBConverter,
)
from modules.file_conversion.converters.editconf_pdb_to_gro import (
    EditconfPDBtoGROConverter,
)
from modules.packmol.solvent_box import PackmolSolventBox
from modules.cache_store.file_cache import FileCache
from typing import List
from modules.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from config.mdp_workflow_config import solvent_workflow
from modules.workflows.base_workflow import BaseWorkflow
from modules.cache_store.solvent_cache import SolventCache

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
solvent_cache = SolventCache()
packmol_solvent_cache = FileCache(name="packmol_solvent_cache")


class SolventEquilibriationWorkflow(BaseWorkflow):
    saved_file_types = ["gro", "top"]
    itp_name = "solvent"

    def __init__(
        self,
        solvent: Solvent,
        parameterised_solvent: GromacsPaths,
        box_size_nm: list,
        temperature: float,
        workflow: FullEquilibrationWorkflow = solvent_workflow,
        output_dir: str = EQUILIBRIATED_SOLVENT_BOX_DIR,
        cleanup: bool = True,
        confirm_temp_dir_deletion: bool = True,
        solvent_cache: SolventCache = solvent_cache,
        packmol_solvent_cache: PickleCache = packmol_solvent_cache,
        verbose: bool = True,
    ):
        self.parameterised_files = parameterised_solvent
        check_directory_exists(TEMP_DIR, make_dirs=True)
        check_directory_exists(LOG_DIR, make_dirs=True)
        self.solvent_cache = solvent_cache
        self.packmol_solvent_cache = packmol_solvent_cache
        self.solvent = solvent
        self.box_size_nm = box_size_nm
        self.temperature = temperature
        self.workflow = workflow
        self.subdir = solvent.name
        self.output_dir = os.path.join(output_dir, self.subdir)
        check_directory_exists(self.output_dir, make_dirs=True)
        self.cleanup = cleanup
        self.verbose = verbose
        self.confirm_temp_dir_deletion = confirm_temp_dir_deletion

    def check_solvent_cache(self, temperature: float) -> Optional[GromacsOutputs]:
        cache_key = self.solvent_cache.get_cache_key(
            solvent=self.solvent, temperature=temperature
        )
        if self.solvent_cache.has_key(cache_key):
            return self.solvent_cache.retrieve_object(cache_key)
        return None

    def _get_packmol_cache_key(self) -> str:
        solvent_name = self.solvent.name
        compressibility = self.solvent.compressibility
        box_size_str = "_".join(map(str, self.box_size_nm))

        return f"{solvent_name}_{compressibility}_{box_size_str}"

    def check_packmol_cache(self) -> Optional[str]:
        cache_key = self._get_packmol_cache_key()
        if self.packmol_solvent_cache.has_key(cache_key):
            return self.packmol_solvent_cache.retrieve_object(cache_key)
        return None

    def _equilibriate(
        self, gro_path: str, top_path: str, temperature: float, itp_path: str
    ) -> Optional[GromacsOutputs]:
        check_file_type(gro_path, "gro")
        check_file_type(top_path, "top")

        outputs = self.check_solvent_cache(temperature)
        if outputs:
            return outputs

        params = [
            {"temp": temperature, "compressibility": self.solvent.compressibility}
        ]
        gro_dir, output_paths = self.workflow.run(
            input_gro_path=gro_path,
            input_topol_path=top_path,
            main_output_dir=self.output_dir,
            temp_output_dir=TEMP_DIR,
            log_dir=LOG_DIR,
            varying_params_list=params,
            files_to_keep=self.saved_file_types,
            subdir=EQUILIBRIATED_OUTPUTS_SUBDIR,
            verbose=self.verbose,
            save_intermediate_edr=True,
            save_intermediate_gro=True,
            save_intermediate_log=True,
        )
        output_paths.itp = itp_path
        cache_key = self.solvent_cache.get_cache_key(
            solvent=self.solvent, temperature=temperature
        )
        solvent_cache.store_object(cache_key, output_paths)
        return output_paths

    def run(self) -> GromacsOutputs:
        parameterised_files = self.parameterised_files
        solvent_box_gro = self.check_packmol_cache()
        if not solvent_box_gro:
            solvent_box_gro = self._create_solvent_box_gro(
                parameterised_files.gro_path,
                PREPROCESSED_PACKMOL_DIR,
                self.box_size_nm,
                self.solvent,
            )
            self.packmol_solvent_cache.store_object(
                self._get_packmol_cache_key(), solvent_box_gro
            )

        reformatted_files = self._process_solvent_files(
            parameterised_files.itp_path,
            solvent_box_gro,
            parameterised_files.top_path,
            new_residue_name=self.solvent.pdb_molecule_name,
            output_itp_dir=self.output_dir,
            output_gro_dir=TEMP_DIR,
            output_topol_dir=TEMP_DIR,
            output_itp_name=self.itp_name,
        )
        output_paths = self._equilibriate(
            gro_path=reformatted_files.gro_path,
            top_path=reformatted_files.top_path,
            temperature=self.temperature,
            itp_path=reformatted_files.itp_path,
        )

        if self.cleanup:
            delete_directory(
                TEMP_DIR, verbose=self.verbose, confirm=self.confirm_temp_dir_deletion
            )
            delete_directory(LOG_DIR, verbose=self.verbose, confirm=False)

        return output_paths

    def _process_solvent_files(
        self,
        input_itp_file: str,
        input_gro_file: str,
        input_top_file: str,
        forcefield: str = "amber99sb-ildn.ff/forcefield.itp",
        new_residue_name: Optional[str] = None,
        output_itp_dir: Optional[str] = None,
        output_gro_dir: Optional[str] = None,
        output_topol_dir: Optional[str] = None,
        output_itp_name: Optional[str] = None,
        output_gro_name: Optional[str] = None,
        output_topol_name: Optional[str] = None,
        parser: GromacsParser = GromacsParser(),
    ) -> GromacsPaths:
        check_file_type(input_itp_file, "itp")
        check_file_type(input_gro_file, "gro")
        check_file_type(input_top_file, "top")

        if new_residue_name:
            if len(new_residue_name) > 5:
                raise ValueError(
                    "Residue name must be 5 characters or less. Gromacs has fixed width for residue names."
                )
        output_itp_file = self._process_solvent_itp(
            input_itp_file=input_itp_file,
            new_residue_name=new_residue_name,
            output_dir=output_itp_dir,
            output_name=output_itp_name,
        )
        output_itp_path = os.path.abspath(output_itp_file)
        gro_handler = get_gro_handler(input_gro_file)
        residue_number = get_residue_number(gro_handler)
        output_top_file = self._prepare_solvent_topol(
            input_top_file=input_top_file,
            residue_number=residue_number,
            forcefield=forcefield,
            new_include_file=output_itp_path,
            new_residue_name=new_residue_name,
            output_name=output_topol_name,
            output_dir=output_topol_dir,
            parser=parser,
            del_posre=True,
            del_defaults=True,
        )
        if new_residue_name:
            gro_handler = rename_residue_name_from_handler(
                gro_handler, new_residue_name
            )

        output_gro_file = prepare_output_file_path(
            input_gro_file, "gro", output_gro_dir, output_gro_name
        )
        output_gro_file = export_gro_handler(gro_handler, output_gro_file, parser)
        paths = GromacsPaths(output_itp_file, output_gro_file, output_top_file)
        return paths

    def _process_solvent_itp(
        self,
        input_itp_file: str,
        new_residue_name: Optional[str] = None,
        output_dir: Optional[str] = None,
        output_name: Optional[str] = None,
        parser: GromacsParser = GromacsParser(),
    ):
        sections = parser.parse(input_itp_file)
        if new_residue_name:
            moleculetype_section = sections["data_moleculetype"]
            moleculetype_section = rename_data_column_content(
                moleculetype_section, "name", new_residue_name
            )
            atoms_section = sections["data_atoms"]
            atoms_section = rename_data_column_content(
                atoms_section, "res", new_residue_name
            )
            sections["data_moleculetype"] = moleculetype_section
            sections["data_atoms"] = atoms_section
        output_itp_path = prepare_output_file_path(
            input_itp_file, "itp", output_dir, output_name
        )
        output_itp_path = parser.export(sections, output_itp_path)
        return output_itp_path

    def _prepare_solvent_topol(
        self,
        input_top_file: str,
        residue_number: int,
        forcefield: str,
        new_include_file: str,
        new_residue_name: Optional[str] = None,  # New single include file for solvent
        output_name: Optional[str] = None,
        output_dir: Optional[str] = None,
        parser: GromacsParser = GromacsParser(),
        del_posre: bool = True,
        del_defaults: bool = True,
    ) -> str:
        residue_number = str(residue_number)
        sections = parser.parse(input_top_file)

        sections = delete_all_include_sections(sections)

        sections["include_solvent_itp"] = create_includes_section(new_include_file)
        sections.move_to_end("include_solvent_itp", last=False)

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
            data_molecules_df = data_molecules_handler.content
            if len(data_molecules_df) != 1:
                raise ValueError("Multiple rows in 'data_molecules' section.")
            if new_residue_name:
                data_molecules_df["Compound"] = new_residue_name
            data_molecules_df["nmols"] = residue_number
            data_molecules_handler.content = data_molecules_df
            data_molecules_section = data_molecules_handler.export()
            sections["data_molecules"] = data_molecules_section

        if del_posre and "conditional_if" in sections:
            del sections["conditional_if"]

        if del_defaults and "data_defaults" in sections:
            del sections["data_defaults"]

        output_path = prepare_output_file_path(
            input_top_file, "top", output_dir, output_name
        )
        output_path = parser.export(sections, output_path)

        return output_path

    def _prepare_solvent_box_name(self, solvent: Solvent, extension: str):
        return f"{solvent.name.lower()}_solvent_box.{extension}"

    def _create_solvent_box_gro(
        self,
        input_gro_file,
        output_dir: str,
        box_size_nm: List[float],
        solvent: Solvent,
        output_name: Optional[str] = None,
        temp_dir: str = TEMP_DIR,
        gro_to_pdb_converter: BaseConverter = EditconfGROtoPDBConverter(),
        pdb_to_gro_converter: BaseConverter = EditconfPDBtoGROConverter(),
        packmol_operation: PackmolSolventBox = PackmolSolventBox(),
    ) -> str:
        output_pdb = gro_to_pdb_converter.run(input_gro_file, temp_dir)
        if not output_name:
            output_name = self._prepare_solvent_box_name(solvent, "gro")

        packmol_output = packmol_operation.run(
            output_pdb,
            output_dir=temp_dir,
            solvent=solvent,
            box_size_nm=box_size_nm,
        )

        output_gro = pdb_to_gro_converter.run(
            packmol_output, output_dir, box_size_nm=box_size_nm, output_name=output_name
        )
        return output_gro
