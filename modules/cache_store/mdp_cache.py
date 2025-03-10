import os
import hashlib
import json
from typing import Dict
from modules.utils.shared.file_utils import (
    check_directory_exists,
    save_content_to_path,
)
import logging

logger = logging.getLogger(__name__)


class MDPCache:
    """
    A class to manage caching of MDP files based on parameters.
    """

    def __init__(self, cache_dir: str):
        """
        Initialize the MDPCache.

        :param cache_dir: Directory to store cached MDP files.
        """
        self.cache_dir = check_directory_exists(cache_dir)
        self.cache_index_path = os.path.join(cache_dir, "cache_index.json")
        self.cache_index = {}
        self._initialize_cache()

    def _initialize_cache(self):
        """
        Load the cache index from file or initialize it.
        """
        if os.path.exists(self.cache_index_path):
            with open(self.cache_index_path, "r") as file:
                self.cache_index = json.load(file)
                logger.info(f"Cache index loaded from {self.cache_index_path}")
        else:
            self.cache_index = {}
            logger.info("Initialized a new cache index.")

    def _save_cache_index(self):
        """
        Save the cache index to file.
        """
        with open(self.cache_index_path, "w") as file:
            json.dump(self.cache_index, file, indent=4)
        logger.info(f"Cache index saved to {self.cache_index_path}")

    def _generate_hash(self, params: Dict[str, str]) -> str:
        """
        Generate a hash key based on the given parameters.

        :param params: Dictionary of MDP parameters.
        :return: A unique hash string for the parameters.
        """
        params_string = json.dumps(params, sort_keys=True)
        hash_key = hashlib.md5(params_string.encode()).hexdigest()
        logger.debug(f"Generated hash key: {hash_key} for params: {params_string}")
        return hash_key

    def _validate_paths(self, template_path: str, output_path: str):
        """
        Validate paths and raise errors for invalid inputs.

        :param template_path: Path to the MDP template file.
        :param output_path: Path to save the generated MDP file.
        """
        if not template_path or not os.path.exists(template_path):
            raise ValueError(f"Invalid template path: {template_path}")

        if not output_path:
            raise ValueError("Output path cannot be None.")

        output_dir = os.path.dirname(output_path)
        if not output_dir:
            raise ValueError(f"Output directory cannot be derived from {output_path}.")

    def _generate_mdp_file(
        self, template_path: str, output_path: str, params: Dict[str, str]
    ):
        """
        Generate an MDP file by replacing placeholders in the template.

        :param template_path: Path to the MDP template file.
        :param output_path: Path to save the generated MDP file.
        :param params: Dictionary of parameters to replace in the template.
        """
        # Validate paths
        self._validate_paths(template_path, output_path)

        # Ensure output directory exists
        check_directory_exists(os.path.dirname(output_path))

        # Load template and replace placeholders
        with open(template_path, "r") as template_file:
            content = template_file.read()

        for key, value in params.items():
            placeholder = f"{{{key}}}"
            content = content.replace(placeholder, str(value))

        # Save the modified content
        save_content_to_path(content.splitlines(keepends=True), output_path)
        logger.info(f"Generated MDP file at {output_path}")

    def get_or_create_mdp(self, template_path: str, params: Dict[str, str]) -> str:
        """
        Retrieve or generate an MDP file based on parameters.

        :param template_path: Path to the MDP template file.
        :param params: Dictionary of parameters for the MDP file.
        :return: Path to the retrieved or newly generated MDP file.
        """
        # Validate cache directory
        check_directory_exists(self.cache_dir)
        logger.info(f"Using cache directory: {self.cache_dir}")
        logger.info(
            f"Checking cache for MDP file with params: {params} with template: {template_path}"
        )

        # Ensure the cache index is initialized
        if not self.cache_index:
            self._initialize_cache()

        # Generate hash key for parameters
        hash_key = self._generate_hash(params)
        mdp_file_path = self.cache_index.get(hash_key)

        if mdp_file_path:
            if os.path.exists(mdp_file_path):
                logger.info(f"MDP file found in cache: {mdp_file_path}")
                return mdp_file_path
            else:
                logger.warning(
                    f"MDP file listed in cache but does not exist: {mdp_file_path}"
                )

        # Derive output path for new MDP file
        mdp_file_path = os.path.join(self.cache_dir, f"{hash_key}.mdp")
        logger.debug(f"Generating new MDP file at: {mdp_file_path}")
        self._generate_mdp_file(template_path, mdp_file_path, params)

        # Update the cache index
        self.cache_index[hash_key] = mdp_file_path
        self._save_cache_index()

        return mdp_file_path
