from abc import ABC, abstractmethod
from typing import Optional, Tuple, List
import subprocess
import logging
from config.constants import LengthUnits

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class CommandLineOperation(ABC):
    def __init__(self):
        pass

    @property
    def step_name(self) -> str:

        return self.__class__.__name__

    @abstractmethod
    def run(self, **kwargs) -> str:
        pass

    def _execute(
        self,
        command: List,
        cwd: Optional[str] = None,
        verbose: bool = False,
        **subprocess_kwargs,
    ) -> str:
        """
        Execute a command using subprocess.run.

        :param command: Command to execute
        :type command: List
        :param cwd: Current working directory, defaults to None
        :type cwd: Optional[str]
        :param verbose: _description_, defaults to False
        :type verbose: bool, optional
        :return: _description_
        :rtype: str
        """
        try:

            self._log_input(command)

            result = subprocess.run(
                command,
                check=True,
                cwd=cwd,
                stdout=subprocess.PIPE if verbose else None,
                stderr=subprocess.PIPE if verbose else None,
                text=True,
                **subprocess_kwargs,
            )

            if verbose:
                self._log_output(result.stdout, result.stderr)

            logger.info(f"Command completed successfully: {self.__class__.__name__}")
        except subprocess.CalledProcessError as e:
            self._handle_error(e, verbose)
            raise

    def _log_input(self, command):
        """
        Log the input command.

        :param command: Command to log
        :type command: List
        """
        logger.info(f"Starting command: {self.__class__.__name__}")
        logger.debug(f"Running command: {' '.join(command)}")

    def _log_output(self, stdout: Optional[str], stderr: Optional[str]):
        """
        Log the output of the command.

        :param stdout: Standard output
        :type stdout: Optional[str]
        :param stderr: Standard error
        :type stderr: Optional[str]
        """
        if stdout:
            logger.debug(f"STDOUT: {stdout}")
        if stderr:
            logger.debug(f"STDERR: {stderr}")

    def _handle_error(self, error: subprocess.CalledProcessError, verbose: bool):
        """
        Handle a CalledProcessError.

        :param error: CalledProcessError to handle
        :type error: subprocess.CalledProcessError
        :param verbose: _description_
        :type verbose: bool
        """
        logger.error(f"Error running  command: {error}")
        if verbose:
            if error.stdout:
                logger.error(f"Command stdout: {error.stdout}")
            if error.stderr:
                logger.error(f"Command stderr: {error.stderr}")
