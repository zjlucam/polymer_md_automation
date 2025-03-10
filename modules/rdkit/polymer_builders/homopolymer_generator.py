from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from rdkit import Chem
from typing import Optional, Dict, List


class HomopolymerGenerator(BasePolymerGenerator):
    def __init__(self):
        super().__init__()
        self.monomer_bead_map: Dict[str, str] = {}  # Map SMILES → Bead Type

    def _generate_polymer_rdkit(self, monomer_smiles: str, num_units: int) -> Chem.Mol:
        """
        Generates a polymer using a flexible list of monomer types and bead assignments.

        :param monomer_smiles: List of monomer SMILES.
        :param num_units: Corresponding number of monomer units per type.
        :return: Final polymer molecule.
        """
        residue_mol, open_sites, bead_type = self._create_monomer_residue(
            monomer_smiles
        )

        polymer, _, idx = self._create_cap_residues(
            residue_mol, open_sites=open_sites, use_open_site=0, add_to_map=True
        )
        # print(self.debug_print_mol(polymer))
        for _ in range(num_units - 2):

            polymer, idx = self._add_monomer_to_polymer(
                polymer,
                monomer=residue_mol,
                prev_end_idx=idx,
                bead_type=bead_type,
                open_sites=open_sites,
            )

        end_bead_type = self._generate_bead_type()
        end_monomer, bead_type, _ = self._create_cap_residues(
            monomer=residue_mol,
            open_sites=open_sites,
            use_open_site=1,
            add_to_map=False,
        )
        polymer, _ = self._add_monomer_to_polymer(
            polymer,
            monomer=end_monomer,
            prev_end_idx=idx,
            bead_type=end_bead_type,
            open_sites=open_sites,
        )

        Chem.SanitizeMol(polymer)
        return polymer

    def _generate_filename(self, monomer_smiles: str, num_units: int) -> str:
        return f"{monomer_smiles.lower()}_{num_units}"

    def generate_polymer(
        self,
        monomer_smiles: str,
        num_units: int,
        output_dir: str,
        output_name: Optional[str] = None,
        uff_optimise: bool = True,
        overwrite: bool = True,
        save: bool = True,
    ) -> str:
        polymer = self._generate_polymer_rdkit(monomer_smiles, num_units)
        polymer = self._finalise_molecule(polymer, uff_optimise=uff_optimise)
        if save:
            if output_name is None:
                output_name = self._generate_filename(monomer_smiles, num_units)
            else:
                output_name = f"{output_name}.pdb"
            output_path = self._save_as_pdb(
                polymer,
                output_dir,
                output_name=output_name,
                overwrite=overwrite,
            )
            return output_path
        else:
            return polymer
