import os
from modules.utils.shared.file_utils import copy_file, rename_file
from config.acpype_config import AcpypeOutputConfig
from config.data_models.output_types import GromacsPaths


def generate_acpype_paths(
    acpype_output_config: AcpypeOutputConfig, directory: str, molecule_name: str
) -> GromacsPaths:
    """
    Generates paths for ACPYPE output files based on the configuration.

    :param acpype_output_config: File configuration for ACPYPE output files.
    :type acpype_output_config: AcpypeOutputConfig
    :param directory: Directory where the files are located.
    :type directory: str
    :param molecule_name: Name of the molecule (ACPYPE molecule name).
    :type molecule_name: str
    :return: Paths to the ACPYPE output files.
    :rtype: GromacsPaths
    """
    return GromacsPaths(
        itp_path=(
            os.path.join(
                directory,
                acpype_output_config.ITP_FILE_NAME_FORMAT.format(
                    molecule_name=molecule_name
                ),
            )
            if acpype_output_config.itp
            else None
        ),
        gro_path=(
            os.path.join(
                directory,
                acpype_output_config.GRO_FILE_NAME_FORMAT.format(
                    molecule_name=molecule_name
                ),
            )
            if acpype_output_config.gro
            else None
        ),
        top_path=(
            os.path.join(
                directory,
                acpype_output_config.TOP_FILE_NAME_FORMAT.format(
                    molecule_name=molecule_name
                ),
            )
            if acpype_output_config.top
            else None
        ),
        posre_path=(
            os.path.join(
                directory,
                acpype_output_config.POSRE_FILE_NAME_FORMAT.format(
                    molecule_name=molecule_name
                ),
            )
            if acpype_output_config.posre
            else None
        ),
    )


def copy_acpype_files(
    acpype_paths: GromacsPaths,
    dest_dir: str,
    delete_original: bool = False,
    skip_if_exists: bool = False,
) -> GromacsPaths:
    """
    Copies files in an GromacsPaths instance to a destination directory.

    :param acpype_paths: Paths to the files to copy
    :type acpype_paths: GromacsPaths
    :param dest_dir: Destination directory to copy the files to
    :type dest_dir: str
    :param delete_original: Flag to delete files or not, defaults to False
    :type delete_original: bool, optional
    :return: Paths to the copied files
    :rtype: GromacsPaths
    """

    copied_files = [
        (
            copy_file(path, dest_dir, delete_original, skip_if_exists=skip_if_exists)
            if path
            else None
        )
        for path in acpype_paths.to_list()
    ]
    return GromacsPaths(*copied_files)


def rename_acpype_paths(
    acpype_paths: GromacsPaths, new_base_name: str, suppress_warning: bool = False
) -> GromacsPaths:
    """
    Renames files in an GromacsPaths instance with a new base name.

    :param acpype_paths: Paths to the files to rename
    :type acpype_paths: GromacsPaths
    :param new_base_name: New base name for the files
    :type new_base_name: str
    :return: Paths to the renamed files
    :rtype: GromacsPaths
    """
    renamed_files = [
        (
            rename_file(path, new_base_name, suppress_warning=suppress_warning)
            if path
            else None
        )
        for path in acpype_paths.to_list()
    ]
    return GromacsPaths(*renamed_files)
