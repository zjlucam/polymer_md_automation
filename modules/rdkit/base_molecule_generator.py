from abc import ABC, abstractmethod
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from modules.utils.shared.file_utils import (
    check_directory_exists,
    check_file_does_not_exist,
)


class BaseMoleculeGenerator(ABC):
    def __init__(self):
        pass

    def _convert_to_rdkit_molecule(self, smiles: str) -> Chem.Mol:
<<<<<<< HEAD
        mol= Chem.MolFromSmiles(smiles)
        AllChem.AddHs(mol)
        return mol
=======
        return Chem.MolFromSmiles(smiles)
>>>>>>> 91758eb (cleaned up)

    def _finalise_molecule(self, mol: Chem.Mol, uff_optimise: bool = True) -> Chem.Mol:
        AllChem.EmbedMolecule(mol)
        if uff_optimise:
            AllChem.UFFOptimizeMolecule(mol)
        AllChem.AddHs(mol)
        return mol

    @abstractmethod
    def _generate_filename(self, **kwargs):
        pass

<<<<<<< HEAD
=======
    def _add_hydrogens(self, molecule: Chem.Mol) -> Chem.Mol:
        return Chem.AddHs(molecule)

>>>>>>> 91758eb (cleaned up)
    def _save_as_pdb(
        self,
        molecule: Chem.Mol,
        output_dir: str,
        output_name: str,
        overwrite: bool = True,
    ) -> str:
        output_basename = output_name + ".pdb"
        output_path = os.path.join(output_dir, output_basename)
        check_directory_exists(output_dir)
        if overwrite:
            suppress_error = True
            delete_file = True
        else:
            suppress_error = False
            delete_file = False
        check_file_does_not_exist(
            output_path, suppress_error=suppress_error, delete_file=delete_file
        )
        with open(output_path, "w") as f:
            f.write(Chem.MolToPDBBlock(molecule))

        self.pdb_path = output_path
        return output_path

    def _get_molar_mass(self, molecule: Chem.Mol, decimal_places: int = 2) -> float:
        molar_mass = AllChem.CalcExactMolWt(molecule)
        return round(molar_mass, decimal_places)
