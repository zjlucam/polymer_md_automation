from modules.gromacs.parsers.handlers.base_handler import (
    BaseHandler,
)
import pandas as pd
import re
from typing import List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class GroHandler(BaseHandler):
    re_pattern = None
    construct_name = "gro"
    suppress = None

    def __init__(
        self,
        expected_columns: List[str] = [
            "Residue Number",
            "Residue Name",
            "Atom Name",
            "Atom Index",
            "X",
            "Y",
            "Z",
            "In-line comments",
        ],
    ):
        super().__init__(store_top_line=True)  # Store the top line
        self.expected_columns = expected_columns
        self.num_atoms = 0
        self.atom_data = []
        self._box_dimensions = None

    @property
    def box_dimensions(self) -> List[float]:
        """
        Returns the box dimensions as a list of three floats.
        """
        return self._box_dimensions

    @box_dimensions.setter
    def box_dimensions(self, box_dimensions: Optional[List[float]]):
        """
        Sets the box dimensions. Validates that the input is a list of three floats.
        Allows setting to None.
        """
        if box_dimensions is None:
            self._box_dimensions = None
            return

        if not isinstance(box_dimensions, list) or len(box_dimensions) != 3:
            raise ValueError("Box dimensions must be a list of three floats.")

        self._box_dimensions = [float(dim) for dim in box_dimensions]

    def process(self, section):
        """
        Process a .gro file section to parse atom lines and extract box dimensions.
        """
        self.section = section
        lines = section.lines  # Do not strip lines to preserve fixed-width format

        # Ensure there is only one top line, concatenating if necessary
        if lines:
            possible_top_lines = []
            while (
                lines and not lines[0].strip().isdigit()
            ):  # Until the atom count is found
                possible_top_lines.append(lines.pop(0).strip())
            self.top_line = (
                " ".join(possible_top_lines)
                if possible_top_lines
                else "default top line"
            )

        else:
            self.top_line = "default top line"

        # Parse number of atoms
        if not lines:
            raise ValueError("Missing number of atoms line in the .gro file.")
        try:
            self.num_atoms = int(lines.pop(0).strip())
        except ValueError:
            raise ValueError("Expected an integer for the number of atoms.")

        # Parse atom data and box dimensions
        for line in lines:
            if len(self.atom_data) < self.num_atoms:
                normalized_tokens = self._normalize_atom_line(line)
                self.atom_data.append(normalized_tokens)
            else:
                # Parse box dimensions
                box_dims = self.parse_box_dimensions(line.strip())
                if self.validate_box_dimensions(box_dims):
                    self.box_dimensions = box_dims
                else:
                    logger.warning(f"Invalid box dimensions: {line.strip()}")
                break

        # Validate that all atoms and box dimensions are parsed correctly
        if len(self.atom_data) != self.num_atoms:
            logger.error(
                f"Expected {self.num_atoms} atoms but parsed {len(self.atom_data)}."
            )

            raise ValueError("Mismatch between expected and parsed number of atoms.")

        if self.box_dimensions is None:
            logger.warning("No valid box dimensions found.")

    def _normalize_atom_line(self, line: str) -> List:
        """
        Normalize a `.gro` atom line using fixed-width parsing.
        Ensures fields are extracted correctly even if misaligned due to spacing.
        """

        comment = None
        if ";" in line:
            content, comment = line.split(";", 1)
            line = content.rstrip()
            comment = comment.strip()

        # Fixed-width parsing
        try:
            residue_number = int(line[0:5].strip())
            residue_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_index = int(line[15:20].strip())
            x = float(line[20:28].strip())
            y = float(line[28:36].strip())
            z = float(line[36:44].strip())
        except ValueError as e:
            raise ValueError(f"Failed to parse line: '{line}'. Details: {e}")

        return [residue_number, residue_name, atom_name, atom_index, x, y, z, comment]

    def _export_content(self) -> List[str]:
        """
        Export the parsed content to `.gro` format lines.
        Ensures adherence to the fixed-width format for both atom lines and box dimensions.
        """
        lines = []

        # Add the atom count as the second line
        actual_num_atoms = len(self.atom_data)

        if self.num_atoms != actual_num_atoms:
            logger.warning(
                f"Mismatch between expected ({self.num_atoms}) and actual ({actual_num_atoms}) atom counts."
            )
            self.num_atoms = actual_num_atoms  # Update to reflect the actual count
        lines.append(f"{self.num_atoms}")

        # Export atom data with strict fixed-width formatting
        for row in self.atom_data:
            content = (
                f"{row[0]:>5}"  # Residue number (right-aligned, width 5)
                f"{row[1]:<5}"  # Residue name (left-aligned, width 5)
                f"{row[2]:>5}"  # Atom name (right-aligned, width 5)
                f"{row[3]:>5}"  # Atom index (right-aligned, width 5)
                f"{row[4]:8.3f}"  # X coordinate (fixed width, 8 chars, 3 decimals)
                f"{row[5]:8.3f}"  # Y coordinate (fixed width, 8 chars, 3 decimals)
                f"{row[6]:8.3f}"  # Z coordinate (fixed width, 8 chars, 3 decimals)
            )
            if row[7]:  # Include in-line comments if present
                content += f" ; {row[7]}"
            lines.append(content)

        # Export box dimensions, ensuring proper formatting
        if self.validate_box_dimensions(self.box_dimensions):
            box_line = f"{self.box_dimensions[0]:10.5f} {self.box_dimensions[1]:10.5f} {self.box_dimensions[2]:10.5f}"
            lines.append(box_line)
        else:
            raise ValueError(
                "Missing or invalid box dimensions; cannot export .gro file."
            )

        return lines

    @staticmethod
    def parse_box_dimensions(line: str) -> Optional[List[float]]:
        """
        Parse a line containing box dimensions.
        Returns:
            List[float]: Parsed box dimensions if valid.
            None: If the line does not contain valid box dimensions.
        """
        try:
            tokens = line.split()
            if len(tokens) == 3 and all(GroHandler._is_float(dim) for dim in tokens):
                return [float(dim) for dim in tokens]
        except ValueError:
            pass
        return None

    @staticmethod
    def validate_box_dimensions(box_dimensions: Optional[List[float]]) -> bool:
        """
        Validate that box dimensions are a list of three positive floats.
        """
        if not box_dimensions or len(box_dimensions) != 3:
            return False
        return all(isinstance(dim, float) and dim > 0 for dim in box_dimensions)

    @staticmethod
    def _is_float(value: str) -> bool:
        """
        Check if a string can be converted to a float.
        """
        try:
            float(value)
            return True
        except ValueError:
            return False

    @property
    def content(self) -> pd.DataFrame:
        """
        Return atom data as a Pandas DataFrame.
        """
        return pd.DataFrame(self.atom_data, columns=self.expected_columns)

    @content.setter
    def content(self, new_content: pd.DataFrame):
        if list(new_content.columns) != self.expected_columns:
            raise ValueError(
                "Columns of the DataFrame do not match the expected format."
            )
        self.atom_data = new_content.values.tolist()
