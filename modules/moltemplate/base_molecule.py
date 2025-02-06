<<<<<<< HEAD
from modules.lammps.parsers.open_mscg_data_parser import OpenMSCGDataParser
from abc import ABC, abstractmethod
from typing import Optional, Dict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BaseMoltemplateMolecule(ABC):

    def __init__(
        self,
        openmscg_parser: OpenMSCGDataParser,
        use_openmscg_mass_map: bool = False,
    ):
        self.premade_mass_map = openmscg_parser.mass_mapping
        self.use_openmscg_mass_map = use_openmscg_mass_map
        self.mass_mapping = {}
        self.lt_block = None

    def _compare_mass_maps(self, replace_using_openmscg_mass_map: bool = False) -> None:
        if replace_using_openmscg_mass_map:
            reference_map = self.premade_mass_map
            compared_map = self.mass_mapping
        else:
            reference_map = self.mass_mapping
            compared_map = self.premade_mass_map

        if not compared_map or not reference_map:
            logger.warning(
                "No mass map provided for comparison. Skipping mass map comparison."
            )
            return

        common_beads = set(reference_map.keys()) & set(compared_map.keys())

        for bead_type in common_beads:
            reference_mass, compared_mass = (
                reference_map[bead_type],
                compared_map[bead_type],
            )

            if round(reference_mass) != round(compared_mass):
                logger.warning(
                    f"Masses for bead type {bead_type} do not match. For mass_map_1: {reference_mass}, for mass_map_2: {compared_mass}"
                )
                if replace_using_openmscg_mass_map:
                    logger.warning(
                        f"Replacing mass for bead type {bead_type} with mass from openmscg mass map."
                    )
                    self.mass_mapping[bead_type] = reference_mass

    @abstractmethod
    def _generate_lt(self) -> str:
        pass

    def generate_lt(self) -> str:
        self.mass_mapping = self._get_mass_mapping()
        self._compare_mass_maps(
            replace_using_openmscg_mass_map=self.use_openmscg_mass_map
        )
        return self._generate_lt()

    @abstractmethod
    def _get_mass_mapping(self) -> Dict[str, float]:
        pass
=======
from modules.lammps.parsers.open_mscg_data_parser import OpenMSCGDataParser
from abc import ABC, abstractmethod
from typing import Optional, Dict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BaseMoltemplateMolecule(ABC):

    def __init__(
        self,
        openmscg_parser: OpenMSCGDataParser,
        use_openmscg_mass_map: bool = False,
    ):
        self.premade_mass_map = openmscg_parser.mass_mapping
        self.use_openmscg_mass_map = use_openmscg_mass_map
        self.mass_mapping = {}
        self.lt_block = None

    def _compare_mass_maps(self, replace_using_openmscg_mass_map: bool = False) -> None:
        if replace_using_openmscg_mass_map:
            reference_map = self.premade_mass_map
            compared_map = self.mass_mapping
        else:
            reference_map = self.mass_mapping
            compared_map = self.premade_mass_map

        if not compared_map or not reference_map:
            logger.warning(
                "No mass map provided for comparison. Skipping mass map comparison."
            )
            return

        common_beads = set(reference_map.keys()) & set(compared_map.keys())

        for bead_type in common_beads:
            reference_mass, compared_mass = (
                reference_map[bead_type],
                compared_map[bead_type],
            )

            if round(reference_mass) != round(compared_mass):
                logger.warning(
                    f"Masses for bead type {bead_type} do not match. For mass_map_1: {reference_mass}, for mass_map_2: {compared_mass}"
                )
                if replace_using_openmscg_mass_map:
                    logger.warning(
                        f"Replacing mass for bead type {bead_type} with mass from openmscg mass map."
                    )
                    self.mass_mapping[bead_type] = reference_mass

    @abstractmethod
    def _generate_lt(self) -> str:
        pass

    def generate_lt(self) -> str:
        self.mass_mapping = self._get_mass_mapping()
        self._compare_mass_maps(
            replace_using_openmscg_mass_map=self.use_openmscg_mass_map
        )
        return self._generate_lt()

    @abstractmethod
    def _get_mass_mapping(self) -> Dict[str, float]:
        pass
>>>>>>> 91758eb (cleaned up)
