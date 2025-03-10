from ruamel.yaml import YAML
from typing import List, Union, Dict, Optional
from modules.utils.shared.file_utils import check_directory_exists
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class OpenMSCGTopolGenerator:
    def __init__(self, map_path: str):
        self.map_data = self._parse_yaml(map_path)
        self.mol_counter = 1
        self.mols_list = []
        self.top_path = None

        self._process_yaml()

    def _parse_yaml(self, map_path: str):
        yaml = YAML()
        with open(map_path, "r") as file:
            data = yaml.load(file)
        return data

    def _get_cgsites(self):
        total_sites = 0

        for entry in self.map_data.get("system"):
            if "repeat" not in entry:
                logger.warning(f"Missing 'repeat' field in system entry: {entry}")
                raise ValueError("Missing 'repeat' field in system entry")
            if "sites" not in entry:
                logger.warning(f"Missing 'repeat' field in system entry: {entry}")
                raise ValueError("Missing 'sites' field in system entry")
            total_sites += entry.get("repeat") * len(entry.get("sites"))

        return total_sites

    def _process_yaml(self) -> None:
        mol_dict = {}
        for entry in self.map_data.get("system"):
            sites_tuple = tuple(tuple(site) for site in entry["sites"])
            if sites_tuple in mol_dict:
                mol_dict[sites_tuple]["count"] += entry["repeat"]
            else:
                mol_dict[sites_tuple] = {
                    "id": self.mol_counter,
                    "count": entry["repeat"],
                    "sites": entry["sites"],
                }
                self.mol_counter += 1

        self.mols_list = [
            {
                "molecule_name": data["id"],
                "sites": data["sites"],
                "repeat_count": data["count"],
            }
            for data in mol_dict.values()
        ]

    def _extract_unique_site_types(
        self, sites: List[List[Union[str, int]]]
    ) -> List[str]:
        return list(dict.fromkeys(site[0] for site in sites))

    def _generate_mol_line(self, sites: List[List[Union[str, int]]]) -> str:
        return f"mol {str(len(sites))} 1"

    def _generate_sitetypes_section(self, sites: List[List[Union[str, int]]]) -> str:
        return "sitetypes\n" + "\n".join(str(site[0]) for site in sites)

    def _get_site_num(self, site_name, sites) -> int:
        return sites[0].index(site_name) + 1

    def _generate_bonds_section(
        self,
        sites: List[List[Union[str, int]]],
    ) -> str:
        bond_list = []
        for i in range(1, len(sites)):
            bond_list.append(f"{i} {i+1}")
        content = f"bonds {str(len(sites) - 1)}"
        if bond_list:
            content += "\n" + "\n".join(bond_list)
        return content

    def _generate_mol_section(
        self, mol_summary: Dict[str, Union[int, List[List[Union[str, int]]]]]
    ):
        mol_line = self._generate_mol_line(sites=mol_summary["sites"])
        sitetypes_lines = self._generate_sitetypes_section(sites=mol_summary["sites"])
        bonds_lines = self._generate_bonds_section(sites=mol_summary["sites"])
        return f"{mol_line}\n{sitetypes_lines}\n{bonds_lines}"

    def _generate_moltypes_section(self):
        moltype_lines = [f"moltypes {str(len(self.mols_list))}"]
        for mol in self.mols_list:
            self._generate_mol_section(mol)
            moltype_lines.append(self._generate_mol_section(mol))
        return "\n".join(moltype_lines)

    def _generate_system_section(self):
        system_lines = [f"system {str(len(self.mols_list))}"]
        for mol in self.mols_list:
            system_lines.append(f"{mol['molecule_name']} {mol['repeat_count']}")
        return "\n".join(system_lines)

    def _generate_cgsites_section(self):
        total_sites = sum(
            len(molecule["sites"] * molecule["repeat_count"])
            for molecule in self.mols_list
        )
        return f"cgsites {str(total_sites)}"

    def _generate_cgtypes_section(self):
        unique_sites = set(
            site[0] for molecule in self.mols_list for site in molecule["sites"]
        )
        return f"cgtypes {str(len(unique_sites))} \n" + "\n".join(unique_sites)

    def _generate_topol(self):
        topol_lines = [
            self._generate_cgsites_section(),
            self._generate_cgtypes_section(),
            self._generate_moltypes_section(),
            self._generate_system_section(),
        ]
        return "\n".join(topol_lines)

    def _save_topol(
        self, content: List[str], output_filename: str, output_dir: Optional[str] = None
    ) -> str:
        if output_dir:
            check_directory_exists(output_dir, make_dirs=True)
            output_filename = f"{output_dir}/{output_filename}"
        with open(output_filename, "w") as file:
            file.write(content)
        logger.info(f"Successfully saved topol file to {output_filename}")
        self.top_path = output_filename
        return output_filename

    def create_topol(self, filename: str, output_dir: Optional[str] = None) -> str:
        filename = f"{filename}.top"
        topol_content = self._generate_topol()
        filepath = self._save_topol(
            content=topol_content, output_filename=filename, output_dir=output_dir
        )
        return filepath
