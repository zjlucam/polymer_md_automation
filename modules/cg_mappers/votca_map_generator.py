import xml.etree.ElementTree as ET
from typing import List, Dict
import xml.etree.ElementTree as ET
from typing import List, Tuple, Dict, Optional
import xml.dom.minidom
import pandas as pd
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.data_handler import DataHandler
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from modules.utils.shared.file_utils import check_directory_exists
import logging
from modules.cg_mappers.base_map_generator import BaseMapGenerator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class VOTCAMapGenerator(BaseMapGenerator):
    """
    Generates a VOTCA-compatible XML mapping file for coarse-grained simulations.
    """

    map_extension = "xml"

    def __init__(
        self,
        polymer: BasePolymerGenerator,
        itp_file_path: Optional[str],
        molecule_name: Optional[str] = None,
        verify_itp: bool = False,
    ):
        """
        Initializes the VOTCA XML generator.

        :param molecule_name: Name of the CG molecule.
        :param bead_mappings: List of dicts containing:
                              - "unique_name": Bead identifier.
                              - "bead_type": Bead type.
                              - "atom_indices": List of atom indices.
                              - "atom_names": Predicted atom names.
        :param bonds: Optional list of bead-pair tuples defining bonds.
        :param angles: Optional list of bead-triples defining angles.
        :param verify_itp: If True, validates against an .itp file.
        :param itp_file_path: Path to the .itp file.
        """

        super().__init__(polymer)
        if not molecule_name:
            self._molecule_name = polymer.res_name
        self.itp_data = None
        self.itp_data = self._parse_itp_data(itp_file_path)

        if verify_itp and itp_file_path:
            self._verify_bonds_against_itp()

    @property
    def molecule_name(self) -> str:
        return self._molecule_name

    @molecule_name.setter
    def molecule_name(self, name: str):
        self._molecule_name = name

    def _parse_itp_data(self, itp_file_path: str) -> pd.DataFrame:
        """
        Parses the .itp file into a DataFrame.
        Assumes `convert_to_itp(itp_file_path)` exists and provides a DataFrame
        with expected headers: `["nr", "type", "resnr", "residue", "atom", ...]`

        :param itp_file_path: Path to the .itp file.
        :return: DataFrame containing residue and atom information.
        """
        parser = GromacsParser()
        sections = parser.parse(itp_file_path)
        atom_section = sections["data_atoms"]
        data_handler = DataHandler()
        data_handler.process(atom_section)
        return data_handler.content

    def _get_atom_metadata(self, atom_index: int) -> Tuple[int, str, str]:
        """
        Retrieves (RESID, RESNAME, ATOMNAME) for a given atom index.

        :param atom_index: Atom index.
        :return: (RESID, RESNAME, ATOMNAME)
        """

        adjusted_index = atom_index + 1
        self.itp_data["nr"] = self.itp_data["nr"].astype(int)
        row = self.itp_data[self.itp_data["nr"] == adjusted_index]
        if row.empty:
            raise ValueError(
                f"[ERROR] Atom index {adjusted_index} not found in .itp file!"
            )

        resname = row["res"].values[0]  # Residue name
        atomname = row["atom"].values[0]  # Atom name

        return 1, resname, atomname  # RESID is always 1ID is always 1

    def _generate_mapping(self, start_index: Optional[int] = None) -> str:
        """
        Saves the mapping data to a VOTCA-compatible XML file.
        """
        root = ET.Element("cg_molecule")
        ET.SubElement(root, "name").text = self.molecule_name
        ET.SubElement(root, "ident").text = self.molecule_name

        topology_elem = ET.SubElement(root, "topology")
        cg_beads_elem = ET.SubElement(topology_elem, "cg_beads")

        # Add bead mappings
        for bead in self.bead_mappings:
            bead_elem = ET.SubElement(cg_beads_elem, "cg_bead")
            ET.SubElement(bead_elem, "name").text = bead["unique_name"]
            ET.SubElement(bead_elem, "type").text = bead["bead_type"]
            ET.SubElement(bead_elem, "mapping").text = (
                f"M{self.bead_mappings.index(bead) + 1}"
            )

            # Corrected: Format bead definitions using RESID:RESNAME:ATOMNAME
            bead_entries = []
            for atom in bead["atom_indices"]:
                resid, resname, atomname = self._get_atom_metadata(atom)
                bead_entries.append(f"{resid}:{resname}:{atomname}")

            ET.SubElement(bead_elem, "beads").text = " ".join(bead_entries)

        # **Reintroduced Bonded Interactions**
        if self.bonds or self.angles:
            cg_bonded_elem = ET.SubElement(topology_elem, "cg_bonded")

            if self.bonds:
                bond_elem = ET.SubElement(cg_bonded_elem, "bond")
                ET.SubElement(bond_elem, "name").text = "bond"
                ET.SubElement(bond_elem, "beads").text = "\n".join(
                    f"{b1} {b2}" for b1, b2 in self.bonds
                )

            if self.angles:
                angle_elem = ET.SubElement(cg_bonded_elem, "angle")
                ET.SubElement(angle_elem, "name").text = "angle"
                ET.SubElement(angle_elem, "beads").text = "\n".join(
                    f"{b1} {b2} {b3}" for b1, b2, b3 in self.angles
                )

        # **Reintroduced Mapping Weights**
        maps_elem = ET.SubElement(root, "maps")
        for i, bead in enumerate(self.bead_mappings):
            map_elem = ET.SubElement(maps_elem, "map")
            ET.SubElement(map_elem, "name").text = f"M{i + 1}"

            # Corrected: Ensure weights are extracted from mapping
            weights = " ".join(map(str, bead["x-weight"]))
            ET.SubElement(map_elem, "weights").text = weights

        xml_str = ET.tostring(root, encoding="utf-8").decode("utf-8")
        formatted_xml = self._prettify_xml(xml_str)

        return formatted_xml

    def _prettify_xml(self, xml_str: str) -> str:
        """Formats XML with proper indentation."""
        dom = xml.dom.minidom.parseString(xml_str)
        return dom.toprettyxml(indent="  ")  # Two-space indentation
