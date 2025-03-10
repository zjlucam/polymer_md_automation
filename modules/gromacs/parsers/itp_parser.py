from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.data_handler import DataHandler
import pandas as pd
from typing import Dict, List, Tuple


class ITPParser:

    def __init__(
        self,
        itp_path: str,
        parser: GromacsParser = GromacsParser(),
        handler=DataHandler(),
    ):
        self.parser = parser
        self.handler = handler
        self.sections = parser.parse(itp_path)

    def _get_section_key(self, section_name: str):
        return f"data_{section_name}"

    def retrieve_handler(self, section_name: str) -> DataHandler:
        section_key = self._get_section_key(section_name)
        section = self.sections[section_key]
        self.handler.process(section)
        return self.handler

    def retrieve_content(self, section_name: str) -> pd.DataFrame:
        retrieved_handler = self.retrieve_handler(section_name)
        return retrieved_handler.content

    def modify_section(self, section_name: str, new_df) -> Dict[str, Dict[str, str]]:
        section_key = self._get_section_key(section_name)
        handler = self.retrieve_handler(section_name)
        handler.content = new_df
        modified_section = handler.export()
        self.sections[section_key] = modified_section
        return self.sections
