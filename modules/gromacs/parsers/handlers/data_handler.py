from modules.gromacs.parsers.handlers.base_handler import (
    BaseHandler,
)
from modules.gromacs.parsers.data_models.section import (
    Section,
)
from modules.gromacs.parsers.registries.section_registry import (
    DataRegistry,
)
import re
import pandas as pd
from typing import List, Dict


class DataHandler(BaseHandler):
    re_pattern = re.compile(r"^\s*\[\s*(.+?)\s*\]\s*$")
    construct_name = "data"
    suppress = None

    def __init__(self, section_registry: DataRegistry = DataRegistry()):
        super().__init__(store_top_line=True)  # Handles the top line ([ section ])
        self.section_registry = section_registry  # Instance of SectionRegistry
        self.expected_headers = []  # Populated during processing
        self.data = []  # Stores raw data rows (with in-line comments)

    def process_line(self, line: str):
        """
        Processes a data line, splitting in-line comments if present.
        """
        if ";" in line:
            content, comment = line.split(";", 1)
            self.data.append(content.strip().split() + [comment.strip()])
        else:
            self.data.append(
                line.split() + [None]
            )  # Add None for missing in-line comment

    def process(self, section: Section):
        """
        Processes the Section object and validates headers using SectionRegistry.
        """
        super().process(section)

        # Use the section.name directly for retrieving headers
        section_name = section.name
        if not section_name:
            raise ValueError("Section name is missing.")

        # Retrieve headers from the registry
        self.expected_headers = self.section_registry.get_headers(section_name)

        # Validate headers
        if self.expected_headers:
            provided_headers = self.expected_headers + ["In-Line Comments"]
            if len(self.data[0]) != len(provided_headers):
                raise ValueError(f"Header validation failed for '{section_name}'.")

    @property
    def content(self) -> pd.DataFrame:
        """
        Returns the data block as a DataFrame.
        """
        if not self.expected_headers:
            raise ValueError("Headers are not defined.")
        columns = self.expected_headers + ["In-Line Comments"]
        return pd.DataFrame(self.data, columns=columns)

    @content.setter
    def content(self, new_content: pd.DataFrame):
        """
        Sets the content using a DataFrame and validates headers.
        """
        if not self.expected_headers:
            raise ValueError("Headers are not defined for this section.")
        self.section_registry.validate_headers(
            self.top_line.strip("[ ]"),
            list(new_content.columns[:-1]),  # Exclude in-line comments column
        )
        self.data = new_content.values.tolist()

    def _export_content(self) -> List[str]:
        """
        Exports the data block as lines, appending in-line comments where applicable.
        """
        lines = []
        for row in self.data:
            content = " ".join(row[:-1])  # All columns except the last
            inline_comment = row[-1]
            if inline_comment:
                lines.append(f"{content} ; {inline_comment}")
            else:
                lines.append(content)
        return lines
