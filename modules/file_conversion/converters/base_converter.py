from abc import ABC, abstractmethod
from modules.command_line_operation import CommandLineOperation
from modules.utils.shared.file_utils import check_file_type
from typing import Optional
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class BaseConverter(CommandLineOperation, ABC):
    input_file_type: str
    output_file_type: str
    program: str

    @property
    def step_name(self) -> str:
        return "conversion"

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if not hasattr(cls, "input_file_type") or not hasattr(cls, "output_file_type"):
            raise TypeError(
                f"Class {cls.__name__} must define 'input_file_type' and 'output_file_type'."
            )

    def __init__(self):
        super().__init__()

    def run(
        self,
        input_file_path: str,
        output_dir: Optional[str] = None,
        verbose: bool = False,
        **kwargs,
    ) -> str:
        logger.info(f"Starting conversion for file: {input_file_path}")
        check_file_type(input_file_path, self.input_file_type)
        return self._run_impl(
            input_file_path=input_file_path,
            output_dir=output_dir,
            verbose=verbose,
            **kwargs,
        )

    @abstractmethod
    def _run_impl(
        self,
        input_file_path: str,
        output_dir: Optional[str] = None,
        verbose: bool = False,
        **kwargs,
    ) -> str:
        pass
