from typing import Optional


handlers = {}


# NOTE: need to figure out what construct type and handler name acc is
class Section:
    def __init__(
        self, construct_name: str, handler_name: str, name: Optional[str] = None
    ):
        self.construct_name = construct_name
        self.handler_name = handler_name
        self.name = name  # The name of the section, extracted from the match
        self.lines = []  # Lines belonging to this section

    def add_line(self, line: str):
        """
        Adds a line to the section's content.

        :param line: The line to add.
        """
        self.lines.append(line)

    def __repr__(self):
        """
        String representation of the Section for debugging purposes.
        """
        return f"Section(construct_type={self.construct_name}, handler={self.handler_name}, name={self.name}, lines={len(self.lines)})"
