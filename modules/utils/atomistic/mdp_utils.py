import os
import shutil
import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def format_temperatures(
    temperatures: List[float], compressibility: str
) -> List[Dict[str, str]]:
    return [
        {"temp": str(temp), "compressibility": str(compressibility)}
        for temp in temperatures
    ]


def generate_file_from_template(
    template_path: str,
    output_path: str,
    replacements: Dict[str, str],
    overwrite: bool = False,
) -> str:
    """
    Generates a file based on a template with multiple key-value replacements.

    :param template_path: Path to the template file to read.
    :type template_path: str
    :param output_path: Path to save the generated file.
    :type output_path: str
    :param replacements: Dictionary of placeholder keys and their replacement values.
    :type replacements: Dict[str, str]
    :param overwrite: Whether to overwrite the output file if it already exists.
    :type overwrite: bool
    :return: The path to the generated file.
    :rtype: str
    """
    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Template file not found: {template_path}")

    if os.path.exists(output_path) and not overwrite:
        logger.info(
            f"Output file already exists and overwrite is disabled: {output_path}"
        )
        return output_path

    with open(template_path, "r") as template_file:
        content = template_file.read()

    for key, value in replacements.items():
        placeholder = f"{{{key}}}"  # Assuming placeholders are in `{key}` format
        if placeholder not in content:
            logger.warning(f"Placeholder '{placeholder}' not found in template.")
        content = content.replace(placeholder, value)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as output_file:
        output_file.write(content)

    logger.info(f"Generated file from template: {output_path}")
    return output_path


def generate_dynamic_filename(
    varying_params: Dict[str, str], extension: Optional[str] = None
) -> str:
    """
    Generate a dynamic filename based on varying parameters.

    :param varying_params: Dictionary of parameter names and values.
    :param extension: File extension to use for the filename.
    :return: A filename string in the format `param1_value1_param2_value2.extension`.
    """
    param_str = "_".join(f"{key}_{value}" for key, value in varying_params.items())
    if not extension:
        return param_str
    else:
        return f"{param_str}.{extension}"
