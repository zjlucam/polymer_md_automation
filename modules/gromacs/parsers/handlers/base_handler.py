from typing import List, Optional
from modules.gromacs.parsers.data_models.section import (
    Section,
)  # Assuming Section is a predefined class
import re


class BaseHandler:
    construct_name: str
    re_pattern: Optional[re.Pattern]
    suppress: Optional[List[str]]

    def __init__(self, store_top_line: bool = False):
        self.store_top_line = (
            store_top_line  # Determines if the top line should be stored
        )
        self.top_line = None  # Stores the first line if store_top_line is True
        self.top_comments = []  # Comments before the main content
        self.bottom_comments = []  # Comments after the main content
        self._content = None  # Stores the main content (data or lines)
        self.section = None

    def __init_subclass__(cls):
        super().__init_subclass__()
        required_attrs = ["construct_name", "re_pattern", "suppress"]
        for attr in required_attrs:
            if not hasattr(cls, attr):
                raise TypeError(f"{cls.__name__} must define '{attr}'")

    @property
    def content(self):
        """
        Gets or sets the content of the handler.
        """
        if self._content is None:
            raise ValueError("Content has not been processed yet.")
        return self._content

    @content.setter
    def content(self, new_content):
        self._content = new_content

    @property
    def first_line(self):
        """
        Gets or sets the first line (if store_top_line is True).
        """
        if not self.store_top_line:
            raise AttributeError("This handler does not support a static top line.")
        return self.top_line

    @first_line.setter
    def first_line(self, new_first_line: str):
        if not self.store_top_line:
            raise AttributeError("This handler does not support a static top line.")
        self.top_line = new_first_line

    def process(self, section: Section):
        """
        Processes a Section object and populates the handler's variables.
        """

        self.section = section
        lines = [line.strip() for line in section.lines if line.strip()]

        # Handle the top line if applicable
        if self.store_top_line and lines:
            self.top_line = lines.pop(0).strip()

        # Process the remaining lines
        for line in lines:
            stripped_line = line.strip()
            if stripped_line.startswith(";"):
                # Add to comments
                if self._content is None:
                    self.top_comments.append(stripped_line)
                else:
                    self.bottom_comments.append(stripped_line)
            else:
                self.process_line(stripped_line)

    def process_line(self, line: str):
        """
        Processes a non-comment line. Should be implemented by child classes.
        """
        raise NotImplementedError(
            "The process_line method must be implemented in child classes."
        )

    def export(self) -> Section:
        """
        Reconstructs and returns a Section object with the processed data.
        """
        lines = []
        if self.store_top_line and self.top_line:
            lines.append(self.top_line)  # Add the top line
        lines.extend(self.top_comments)  # Add top comments
        lines.extend(self._export_content())  # Add content
        lines.extend(self.bottom_comments)  # Add bottom comments
        self.section.lines = lines
        return self.section

    def _export_content(self) -> List[str]:
        """
        Exports the main content as lines. Must be implemented by child classes.
        """
        raise NotImplementedError(
            "The _export_content method must be implemented in child classes."
        )
