from modules.gromacs.parsers.handlers.base_handler import (
    BaseHandler,
)
import pandas as pd
from typing import List


class DefaultHandler(BaseHandler):
    re_pattern = None
    construct_name = "default"
    suppress = None

    def __init__(self):
        super().__init__(store_top_line=False)  # No static top line for #include
        self.lines = None  # Stores the #include line

    def process_line(self, line: str):
        self.lines = line.strip()

    @property
    def content(self) -> str:
        return self.lines

    @content.setter
    def content(self, new_content: str):
        self.lines = new_content

    def _export_content(self) -> List[str]:
        return self.lines
