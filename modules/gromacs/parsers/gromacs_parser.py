from collections import OrderedDict
from modules.gromacs.parsers.data_models.section import (
    Section,
)
from modules.gromacs.parsers.handlers.base_handler import (
    BaseHandler,
)
from modules.gromacs.parsers.handlers.gro_handler import (
    GroHandler,
)
from modules.gromacs.parsers.handlers.default_handler import (
    DefaultHandler,
)
from typing import Dict, List, Optional, Tuple, OrderedDict
from modules.gromacs.parsers.registries.handler_registry import (
    HandlerRegistry,
    handler_registry,
)


class GromacsParser:
    def __init__(self, handler_registry: HandlerRegistry = handler_registry):
        self.handler_registry = handler_registry
        self.suppressed_constructs: Optional[List[str]] = None

    # NOTE: make this more robust

    def parse(self, filepath: str) -> OrderedDict[str, Section]:
        sections: OrderedDict[str, Section] = OrderedDict()
        current_section = Section(construct_name=None, handler_name=None)

        with open(filepath, "r") as file:
            for line in file:
                line = line.rstrip("\n")

                construct_name, handler_name, name = self._match_line(line)

                if construct_name:
                    is_data_handler = construct_name == "data"
                    key = self._generate_key(
                        sections,
                        current_section.construct_name,
                        current_section.name,
                    )
                    sections[key] = current_section

                    current_section = Section(
                        construct_name=construct_name,
                        handler_name=handler_name,
                        name=name,
                    )

                current_section.add_line(line)

        if current_section.lines:
            key = self._generate_key(
                sections,
                current_section.construct_name,
                current_section.name,
            )
            sections[key] = current_section

        return sections

    def _generate_key(
        self,
        sections: OrderedDict[str, Section],
        construct_type: str,
        name: Optional[str],
    ) -> str:
        """
        Generate a unique key for a section or handler. If it's a DataHandler,
        append the name to the key; otherwise, use only the handler name.

        Args:
            sections (OrderedDict[str, Section]): Current sections in the parser.
            construct_type (str): Type of the construct (e.g., 'gro_file').
            name (Optional[str]): Name of the section (if applicable).
            handler_type (str): The type of the handler to determine if it's a DataHandler.

        Returns:
            str: Unique key for the section or handler.
        """
        # Identify if it's a DataHandler based on handler type
        is_data_handler = construct_type == "data"

        # Append the name only for DataHandler constructs
        if is_data_handler:
            base_key = f"{construct_type}_{name or 'no_name'}"
        else:
            base_key = construct_type

        if base_key in sections:
            # Append an incremented suffix to ensure uniqueness
            index = 1
            while f"{base_key}_{index}" in sections:
                index += 1
            return f"{base_key}_{index}"

        return base_key

    @property
    def available_handlers(self) -> Dict[str, BaseHandler]:
        """
        Return handlers from the registry, excluding those in suppressed constructs.
        """
        if not self.suppressed_constructs:
            return self.handler_registry._handlers
        return {
            handler_name: handler_class
            for handler_name, handler_class in self.handler_registry._handlers.items()
            if handler_name not in self.suppressed_constructs
        }

    def _match_line(self, line: str) -> Tuple[str, str, str]:
        """
        Matches a line to a construct and returns its type, name, and handler name.
        Filters out suppressed constructs.
        """
        for handler_name, handler_class in self.available_handlers.items():
            if handler_class.re_pattern is None:
                continue

            match = handler_class.re_pattern.match(line)
            if match:
                name = match.group(1) if match.groups() else None
                self.suppressed_constructs = handler_class.suppress
                return handler_class.construct_name, handler_name, name

        # Fallback to DefaultHandler
        return None, DefaultHandler.construct_name, None

    def export(self, sections: OrderedDict[str, Section], output_filepath: str) -> str:
        """
        Export all sections back to a `.gro` file.

        Args:
            sections (OrderedDict[str, Section]): Updated sections to export.
            output_filepath (str): Path to the output `.gro` file.
        """
        with open(output_filepath, "w") as file:
            for section in sections.values():
                file.write("\n".join(section.lines) + "\n")

        return output_filepath
