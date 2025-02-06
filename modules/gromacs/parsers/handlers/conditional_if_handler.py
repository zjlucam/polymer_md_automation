from modules.gromacs.parsers.handlers.base_handler import (
    BaseHandler,
)
from typing import List
import re


class ConditionalIfHandler(BaseHandler):
    construct_name = "conditional_if"
    re_pattern = re.compile(r"^\s*#\s*(ifdef|ifndef)\s+.*$")
    suppress = ["include"]

    def __init__(self):
        super().__init__(
            store_top_line=True
        )  # The top line stores the conditional start (e.g., #ifdef)
        self.conditional_content = []  # Lines within the conditional block

    def process_line(self, line: str):
        """
        Processes lines within the conditional block.
        """
        self.conditional_content.append(line)

    @property
    def content(self) -> List[str]:
        """
        Returns the lines within the conditional block.
        """
        return self.conditional_content

    @content.setter
    def content(self, new_content: List[str]):
        """
        Sets the content of the conditional block.
        """
        self.conditional_content = new_content

    def _export_content(self) -> List[str]:
        """
        Exports the lines within the conditional block, surrounded by #ifdef/#endif.
        """
        return self.conditional_content
