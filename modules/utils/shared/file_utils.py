import os
import shutil
import logging
from typing import List, Callable, TypeVar, Any

# from typing import List, Callable, TypeVar, Any, ParamSpec
from pathlib import Path
from functools import wraps
from typing import Optional
import logging

logger = logging.getLogger(__name__)


def check_file_exists(file_path: str, suppress_error: bool = False) -> Optional[str]:
    """
    Checks that a file exists in the correct location

    :param file_path: Relative or absolute path to the file
    :type file_path: str
    :param suppress_error: If True, suppresses the FileNotFound error when file doesn't exist, and logs a warning instead, function will return None, defaults to False
    :type suppress_error: bool
    :raises FileNotFoundError: If the file is not found
    :return: The file path, if found, else None
    :rtype: Optional[str]
    """

    if not os.path.exists(file_path):
        message = f"File not found: {file_path}"
        if suppress_error:
            logger.warning(message)
            return None
        logger.error(message)
        raise FileNotFoundError(message)

    logger.info(f"File found: {file_path}")
    return file_path


def check_directory_exists(
    directory_path: str, make_dirs: bool = True
) -> Optional[str]:
    """
    Checks that a directory exists in the correct location

    :param directory_path: Relative or absolute path to the directory
    :type directory_path: str
    :param make_dirs: If True, missing directory is created, defaults to True
    :type make_dirs: bool, optional
    :raises FileNotFoundError: If the directory is not found
    :return: The directory path, if found, else None
    :rtype: Optional[str]
    """
    if not os.path.exists(directory_path):
        if make_dirs:
            logger.info(f"Creating directory: {directory_path}")

            os.makedirs(directory_path, exist_ok=True)

            logger.info(f"Directory created: {directory_path}")
            return directory_path
        else:
            logger.error(f"Directory not found: {directory_path}")

            raise FileNotFoundError(f"Directory not found: {directory_path}")
    else:
        logger.info(f"Directory found: {directory_path}")

        return directory_path


def move_files(
    file_paths: List[str], target_directory: str, overwrite: bool = False
) -> List[str]:
    os.makedirs(target_directory, exist_ok=True)
    moved_files = []
    for file_path in file_paths:
        if not os.path.isfile(file_path):
            logger.warning(f"Warning: File not found - {file_path}")
            continue

        file_name = os.path.basename(file_path)
        destination = os.path.join(target_directory, file_name)

        if os.path.exists(destination) and not overwrite:
            logger.info(f"Skipping: {destination} already exists (overwrite=False)")
            moved_files.append(destination)
            continue

        try:
            shutil.move(file_path, destination)
            moved_files.append(destination)
            logger.info(f"Moved: {file_path} -> {destination}")
        except Exception as e:
            logger.error(f"Error moving {file_path}: {e}")
    return moved_files


def get_file_contents(file_path: str, get_contents: bool = True) -> Optional[List[str]]:
    """
    Checks that the contents of a file is not empty, and retrieves contents if :param get_contents: is True

    :param file_path: Relative or absolute path to the file
    :type file_path: str
    :param get_contents: If True, returns the contents of the file as a list of strings, defaults to False
    :type get_contents: bool, optional
    :raises ValueError: If the file is empty
    :return: The contents of the file, if found and :param get_contents: is True, else None
    :rtype: Optional[List[str]]
    """
    with open(file_path, "r") as file:
        lines = file.readlines()

    if not lines:
        logger.error(f"File is empty: {file_path}")
        raise ValueError(f"File is empty: {file_path}")

    logger.info(f"File content checked: {file_path} is not empty")
    if not get_contents:
        return None
    else:
        return lines


def save_content_to_directory(
    content: List[str],
    file_name: str,
    file_extension: str,
    output_dir: str,
    suppress_error: bool = True,
) -> str:
    """
    Saves content to a file

    :param content: Content to save
    :type content: List[str]
    :param file_name: Name of the file without the extension
    :type file_name: str
    :param file_extension: Name of file extension, accepted format is "pdb", not ".pdb"
    :type file_extension: str
    :param output_dir: Directory to save the file
    :type output_dir: str
    :param suppress_error: suppresses the error when the file already exists, and logs a warning instead, function will return None, defaults to True
    :type suppress_error: bool, optional
    :return: The path to the saved file
    :rtype: str
    """

    output_dir = check_directory_exists(output_dir)
    output_file_path = os.path.join(output_dir, f"{file_name}.{file_extension}")
    output_file_path = check_file_does_not_exist(
        output_file_path, suppress_error=suppress_error
    )
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
    with open(output_file_path, "w") as file:
        file.writelines(content)
    logger.info(f"Saved content to {output_file_path}")
    return output_file_path


def save_content_to_path(
    content: List[str],
    output_path: str,
    suppress_error: bool = True,
) -> str:
    """
    Saves content to a file

    :param content: Content to save
    :type content: List[str]
    :param output_path: Path to save the file
    :type output_path: str
    :param suppress_error: suppresses the error when the file already exists, and logs a warning instead, function will return None, defaults to True
    :type suppress_error: bool, optional
    :return: The path to the saved file
    :rtype: str
    """

    check_file_does_not_exist(output_path, suppress_error=suppress_error)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as file:
        file.writelines(content)
    logger.info(f"Saved content to {output_path}")
    return output_path


def check_file_type(file_path: str, expected_file_type: str) -> None:
    """
    Validates that the file types of the input match the expected input_file_type

    :param file_path: Relative or absolute path to the file
    :type file_path: str
    :param expected_file_type: The expected file type
    :type expected_file_type: str
    :raises ValueError: If the file type does not match the expected file type
    """
    observed_file_type = os.path.splitext(file_path)[1].lstrip(".")

    if observed_file_type != expected_file_type:
        message = (
            f"Validation failed: Expected input file of type "
            f"'{expected_file_type}', but got '{observed_file_type}'."
        )
        logger.error(message)
        raise ValueError(message)

    logger.info(f"Validation passed: Input file is of type '.{expected_file_type}'")


def batch_copy_file(
    files_to_move: List[Optional[str]],
    dest_dir: str,
    delete_original: bool = False,
) -> List[Optional[str]]:
    """
    Copies specified files to a destination directory. Optionally deletes the originals.

    :param files_to_move: List of specific file names to copy, file name can be None, but will return None
    :type files_to_move: List[Optional[str]]
    :param dest_dir: The directory to copy files to
    :type dest_dir: str
    :param delete_original: If True, deletes the original files after copying, defaults to False
    :type delete_original: bool, optional
    :return: List of file paths in the destination directory
    :rtype: List[Optional[str]]
    """

    copied_files = []

    for file in files_to_move:
        copied_file = copy_file(
            file_path=file, dest_dir=dest_dir, delete_original=delete_original
        )
        copied_files.append(copied_file)

    return copied_files


def copy_and_rename(
    file_path: Optional[str],
    dest_dir: str,
    new_name: Optional[str] = None,
    delete_original: bool = False,
    replace_if_exists: bool = True,
) -> Optional[str]:
    """
    Copies a file to a destination directory and optionally renames it. Optionally deletes the original.
    Optionally skips copying if the file already exists.

    :param file_path: Relative or absolute path to the file to copy, can be None (will return None).
    :type file_path: Optional[str]
    :param dest_dir: The directory to copy the file to.
    :type dest_dir: str
    :param new_name: New name for the copied file (without extension). If None, keeps the original name.
    :type new_name: Optional[str]
    :param delete_original: If True, deletes the original file after copying, defaults to False.
    :type delete_original: bool, optional
    :param skip_if_exists: If True, skips copying if the file already exists in the destination directory.
    :type skip_if_exists: bool, optional
    :return: The path to the copied (and optionally renamed) file, if successful, else None.
    :rtype: Optional[str]
    """

    if file_path is None:
        logger.info("File is None, skipping.")
        return None

    check_file_exists(file_path)
    check_directory_exists(dest_dir)

    # Determine the new file name
    original_extension = os.path.splitext(file_path)[1]
    new_file_name = (
        f"{new_name}{original_extension}" if new_name else os.path.basename(file_path)
    )
    dest_file_path = os.path.join(dest_dir, new_file_name)

    if os.path.exists(dest_file_path):
        if not replace_if_exists:
            logger.info(f"File already exists at {dest_file_path}. Skipping copy.")
            return dest_file_path
        else:
            logger.info(
                f"File already exists at {dest_file_path}. Deleting existing file."
            )
            os.remove(dest_file_path)

    # Perform the copy operation
    dest_file = shutil.copy2(file_path, dest_file_path)
    logger.info(f"Copied {file_path} to {dest_file}.")

    if delete_original:
        os.remove(file_path)
        logger.info(f"Deleted original file: {file_path}")

    return dest_file


def copy_file(
    file_path: Optional[str],
    dest_dir: str,
    delete_original: bool = False,
    skip_if_exists: bool = True,
) -> Optional[str]:
    """
    Copies a file to a destination directory. Optionally deletes the original. Optionally skips copying if the file already exists.

    :param file_path: Relative or absolute path to the file to copy, can be None (will return None).
    :type file_path: Optional[str]
    :param dest_dir: The directory to copy the file to.
    :type dest_dir: str
    :param delete_original: If True, deletes the original file after copying, defaults to False.
    :type delete_original: bool, optional
    :param skip_if_exists: If True, skips copying if the file already exists in the destination directory.
    :type skip_if_exists: bool, optional
    :return: The path to the copied file, if successful, else None.
    :rtype: Optional[str]
    """

    if file_path is None:
        logger.info("File is None, skipping.")
        return None

    check_file_exists(file_path)
    check_directory_exists(dest_dir)

    dest_file_path = os.path.join(dest_dir, os.path.basename(file_path))

    if skip_if_exists and os.path.exists(dest_file_path):
        logger.info(f"File already exists at {dest_file_path}. Skipping copy.")
        return dest_file_path

    # Perform the copy operation
    dest_file = shutil.copy2(file_path, dest_dir)
    logger.info(f"Copied {file_path} to {dest_file}.")

    if delete_original:
        os.remove(file_path)
        logger.info(f"Deleted original file: {file_path}")

    return dest_file


def batch_rename_to_same(file_paths: List[str], new_name: str) -> List[Optional[str]]:
    """
    Renames multiple files to have the same name, preserving directories and extensions.

    :param file_paths: List of specific file names to rename, file name can be None, but will return None
    :type file_paths: List[str]
    :param new_name: The new name for ALL the files, if you want to specify names, please use batch_rename_to_list()
    :type new_name: str
    :return: List of renamed file paths
    :rtype: List[Optional[str]]
    """
    renamed_files = []

    for file in file_paths:
        renamed_file = rename_file(file_path=file, new_name=new_name)
        renamed_files.append(renamed_file)

    return renamed_files


def batch_rename_to_list(
    file_paths: List[str], new_names: List[str]
) -> List[Optional[str]]:
    """
    Renames multiple files to have the names specified in new_names, preserving directories and extensions.

    :param file_paths: List of specific file names to rename, file name can be None, but will return None
    :type file_paths: List[str]
    :param new_names: List of new names for the files, if you want to rename it all to the same name, please use batch_rename_to_same()
    :type new_names: List[str]
    :return: List of renamed file paths
    :rtype: List[Optional[str]]
    """

    renamed_files = []

    for new_name, file in zip(new_names, file_paths):
        renamed_file = rename_file(file_path=file, new_name=new_name)
        renamed_files.append(renamed_file)

    return renamed_files


def rename_file(
    file_path: Optional[str], new_name: str, suppress_warning: bool = False
) -> Optional[str]:
    """
    Renames a file to a new name.

    :param file_path: Relative or absolute path to the file to rename, can be None (will return None)
    :type file_path: Optional[str]
    :param new_name: The new name for the file
    :type new_name: str
    :param suppress_warning: If True, suppresses the warning when the file already exists, defaults to False
    :type suppress_warning: bool, optional
    :return: The new file path
    :rtype: str
    """
    if file_path is None:
        logger.info("File is None, skipping.")
        return None

    file_path = check_file_exists(file_path)
    path = Path(file_path)

    new_name_with_file_extension = new_name + path.suffix
    new_file_path = path.with_name(new_name_with_file_extension)

    check_file_does_not_exist(new_file_path, suppress_error=suppress_warning)
    os.rename(file_path, new_file_path)

    logger.info(f"Renamed {file_path} to {new_file_path}")
    return str(new_file_path)


def determine_file_name(
    file_path: str,
    new_file_name: Optional[str] = None,
    new_file_extension: Optional[str] = None,
):
    """
    Determines the new file name based on the input file path, with optional new file name and file extension.

    :param file_path: Relative or absolute path to the file
    :type file_path: str
    :param new_file_name: New file name, if none, it will use original, defaults to None
    :type new_file_name: Optional[str], optional
    :param new_file_extension: New file extension, if output file type is different, defaults to None
    :type new_file_extension: Optional[str], optional
    :return: New file name
    :rtype: str
    """
    path = Path(file_path)

    # Use the provided new file extension or fallback to the original one
    file_extension = f".{new_file_extension}" if new_file_extension else path.suffix

    # Use the provided new file name or fallback to the original stem (name without extension)
    file_name = new_file_name if new_file_name else path.stem

    return file_name + file_extension


def check_file_does_not_exist(
    file_path: str = None, suppress_error: bool = False, delete_file: bool = False
) -> None:
    """
    Checks that a file does not exist in the correct location

    :param file_path: Relative or absolute path to the file
    :type file_path: str
    :param suppress_error: If True, suppresses the error when the file already exists, defaults to False
    :type suppress_error: bool, optional
    :raises FileExistsError: If the file already exists
    """
    if file_path is None:
        logger.info("File is None, skipping.")
        return

    if os.path.exists(file_path):
        if suppress_error:
            logger.warning(f"A file with the name '{file_path}' already exists.")
            if delete_file:
                os.remove(file_path)
                logger.info(f"Deleted existing file: {file_path}")
            return
        raise FileExistsError(f"A file with the name '{file_path}' already exists.")

    logger.info(f"File does not exist: {file_path}")


def delete_directory(directory_path: str, verbose: bool = False, confirm: bool = False):
    """
    Deletes the specified directory and its contents.

    :param directory_path: Path to the directory to delete.
    :param verbose: If True, prints detailed information about the deletion process.
    :param confirm: If True, prompts the user for confirmation before deleting.
    :raises FileNotFoundError: If the directory does not exist.
    :raises ValueError: If the path provided is not a directory.
    """
    if not os.path.exists(directory_path):
        raise FileNotFoundError(f"The directory '{directory_path}' does not exist.")

    if not os.path.isdir(directory_path):
        raise ValueError(f"The path '{directory_path}' is not a directory.")

    if confirm:
        user_input = (
            input(f"Are you sure you want to delete '{directory_path}'? (y/n): ")
            .strip()
            .lower()
        )
        if user_input not in ["y", "yes"]:
            if verbose:
                logger.info(f"Deletion of '{directory_path}' canceled by user.")
            return

    try:
        shutil.rmtree(directory_path)
        if verbose:
            logger.info(f"Directory '{directory_path}' has been deleted successfully.")
    except Exception as e:
        raise RuntimeError(f"Failed to delete directory '{directory_path}': {e}")


# Define type variables for the input parameters and return type of the wrapped function
# P = ParamSpec("P")  # Represents the parameters of the wrapped function
# R = TypeVar("R")  # Represents the return type of the wrapped function


def directory_exists_check_wrapper(dir_arg_index: int, make_dirs: bool = True):
    # -> Callable[[Callable[P, R]], Callable[P, R]]
    """
    A wrapper to ensure a directory exists, optionally creating it if missing.

    :param dir_arg_index: The index of the directory path argument in the wrapped function's arguments.
    :type dir_arg_index: int
    :param make_dirs: If True, create the directory if it does not exist. Defaults to True.
    :type make_dirs: bool
    :return: A decorator that wraps the provided function.
    :rtype: Callable[[Callable], Callable]
    """

    # def decorator(func: Callable[P, R]) -> Callable[P, R]:
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            dir_path = (
                args[dir_arg_index]
                if len(args) > dir_arg_index
                else kwargs.get("directory_path")
            )
            if dir_path is None:
                logger.error("Directory path argument not provided.")
                raise ValueError("Directory path argument not provided.")

            if not os.path.exists(dir_path):
                if make_dirs:
                    logger.info(f"Creating directory: {dir_path}")
                    os.makedirs(dir_path, exist_ok=True)
                    logger.info(f"Directory created: {dir_path}")
                else:
                    logger.error(f"Directory not found: {dir_path}")
                    raise FileNotFoundError(f"Directory not found: {dir_path}")
            else:
                logger.info(f"Directory found: {dir_path}")

            return func(*args, **kwargs)

        return wrapper

    return decorator


def file_does_not_exist_check_wrapper(
    file_arg_index: int, suppress_error: bool = False
):
    # -> Callable[[Callable[P, R]], Callable[P, R | None]]
    """
    A wrapper to ensure a file does not exist before calling the wrapped function.

    :param file_arg_index: The index of the file path argument in the wrapped function's arguments.
    :type file_arg_index: int
    :param suppress_error: If True, suppress the error when the file already exists. Defaults to False.
    :type suppress_error: bool
    :return: A decorator that wraps the provided function.
    :rtype: Callable[[Callable], Callable]
    """

    # def decorator(func: Callable[P, R]) -> Callable[P, R | None]:
    def decorator(func):
        @wraps(func)
        # def wrapper(*args: P.args, **kwargs: P.kwargs) -> R | None:
        def wrapper(*args, **kwargs):
            file_path = (
                args[file_arg_index]
                if len(args) > file_arg_index
                else kwargs.get("file_path")
            )
            if file_path is None:
                logger.info("File path argument not provided, skipping.")
                return func(*args, **kwargs)

            if os.path.exists(file_path):
                if suppress_error:
                    logger.warning(f"File already exists: {file_path}")
                    return None
                else:
                    logger.error(f"File already exists: {file_path}")
                    raise FileExistsError(f"File already exists: {file_path}")

            logger.info(f"File does not exist: {file_path}")
            return func(*args, **kwargs)

        return wrapper

    return decorator


def file_exists_check_wrapper(file_arg_index: int, suppress_error: bool = False):
    # -> Callable[[Callable[P, R]], Callable[P, R | None]]
    """
    A wrapper to check if a file exists before calling the wrapped function.

    :param file_arg_index: The index of the file path argument in the wrapped function's arguments.
    :type file_arg_index: int
    :param suppress_error: If True, suppresses the FileNotFound error when the file doesn't exist,
                           logs a warning instead, and skips calling the wrapped function. Defaults to False.
    :type suppress_error: bool
    :return: A decorator that wraps the provided function.
    :rtype: Callable[[Callable], Callable]
    """

    # def decorator(func: Callable[P, R]) -> Callable[P, R | None]:
    def decorator(func):
        @wraps(func)
        # def wrapper(*args: P.args, **kwargs: P.kwargs) -> R | None:
        def wrapper(*args, **kwargs):
            file_path = (
                args[file_arg_index]
                if len(args) > file_arg_index
                else kwargs.get("file_path")
            )
            if file_path is None:
                logger.error("File path argument not provided.")
                raise ValueError("File path argument not provided.")

            if not os.path.exists(file_path):
                message = f"File not found: {file_path}"
                if suppress_error:
                    logger.warning(message)
                    return None
                logger.error(message)
                raise FileNotFoundError(message)

            logger.info(f"File found: {file_path}")
            return func(*args, **kwargs)

        return wrapper

    return decorator


def file_type_check_wrapper(file_arg_index: int, expected_file_type: str):
    # -> Callable[[Callable[P, R]], Callable[P, R]]
    """
    A wrapper to validate the file type of a specific argument before calling the function.

    :param file_arg_index: The index of the file path argument in the wrapped function's arguments.
    :type file_arg_index: int
    :param expected_file_type: The expected file type (e.g., 'mol2', 'pdb').
    :type expected_file_type: str
    :return: A decorator that wraps the provided function.
    :rtype: Callable[[Callable], Callable]
    """

    # def decorator(func: Callable[P, R]) -> Callable[P, R]:
    def decorator(func):
        @wraps(func)
        # def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        def wrapper(*args, **kwargs):
            file_path = (
                args[file_arg_index]
                if len(args) > file_arg_index
                else kwargs.get("file_path")
            )
            if file_path is None:
                logger.error("File path argument not provided.")
                raise ValueError("File path argument not provided.")

            # Check file extension
            observed_file_type = os.path.splitext(file_path)[1].lstrip(".")
            if observed_file_type != expected_file_type:
                message = (
                    f"Validation failed: Expected input file of type "
                    f"'{expected_file_type}', but got '{observed_file_type}'."
                )
                logger.error(message)
                raise ValueError(message)

            logger.info(
                f"Validation passed: Input file is of type '.{expected_file_type}'"
            )

            return func(*args, **kwargs)

        return wrapper

    return decorator


@file_exists_check_wrapper(file_arg_index=0)
def construct_output_file_path(
    file_path: str,
    new_output_dir: Optional[str] = None,
    new_file_name: Optional[str] = None,
    new_file_extension: Optional[str] = None,
    suppress_warning: bool = False,
    expected_file_type: Optional[str] = None,
) -> str:
    """
    Constructs a new file path based on the input file path, with optional new output directory, file name, and file extension.

    :param file_path: Relative or absolute path to the file
    :type file_path: str
    :param new_output_dir: Output directory, if not specified, will use original output file directory, defaults to None
    :type new_output_dir: Optional[str], optional
    :param new_file_name: New file name, if none, it will use original, defaults to None
    :type new_file_name: Optional[str], optional
    :param new_file_extension: New file extension, if output file type is different, defaults to None
    :type new_file_extension: Optional[str], optional
    :param suppress_warning: If True, suppresses the warning when the file already exists, defaults to False
    :type suppress_warning: bool, optional
    :param expected_file_type: Optional file type check, defaults to None
    :type expected_file_type: Optional[str], optional
    :return: Output file path
    :rtype: str
    """

    path = Path(file_path)

    if expected_file_type:
        check_file_type(file_path, expected_file_type)

    if new_output_dir:
        output_dir = new_output_dir
    else:
        output_dir = os.path.dirname(file_path)
    check_directory_exists(output_dir)

    output_file_name = determine_file_name(file_path, new_file_name, new_file_extension)
    output_file_path = os.path.join(output_dir, output_file_name)

    check_file_does_not_exist(output_file_path, suppress_error=suppress_warning)
    return output_file_path


def overwrite_directory(directory_path: str) -> None:
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)
    os.makedirs(directory_path, exist_ok=True)


import os
from typing import Dict


def generate_file_from_template(
    template_path: str,
    output_path: str,
    replacements: Dict[str, str],
    overwrite: bool = False,
) -> None:
    """
    Reads a template file, replaces placeholders with provided values,
    and writes the modified content to an output file.

    :param template_path: Path to the template file to read.
    :type template_path: str
    :param output_path: Path to the output file to write.
    :type output_path: str
    :param replacements: A dictionary of placeholder strings to their replacement values ["placeholder" : "replacement].
    :type replacements: Dict[str, str]
    :param overwrite: If True, overwrite the output file if it already exists. Defaults to False.
    :type overwrite: bool, optional
    """
    file_exists = check_file_exists(output_path, suppress_error=True)

    if file_exists and not overwrite:
        logger.info(f"Output file already exists at {output_path}, skipping creation.")
        return

    with open(template_path, "r") as file:
        content = file.read()

    for placeholder, value in replacements.items():
        content = content.replace(f"{{{placeholder}}}", str(value))

    with open(output_path, "w") as file:
        file.write(content)

    logger.info(f"Created file from template at {output_path}")


def prepare_output_file_path(
    input_file: str,
    output_extension: str,
    output_dir: Optional[str] = None,
    output_name: Optional[str] = None,
) -> str:
    """
    Constructs an output file path based on the input file path, output directory, and new output name.

    :param input_file: Path to the input file.
    :param output_dir: Directory to save the output file. If None, uses input file directory.
    :param output_name: New file name without extension. If None, uses input file name.
    :return: Full path to the output file with the specified extension.
    """
    input_path = Path(input_file)
    output_dir = Path(output_dir) if output_dir else input_path.parent
    output_name = output_name if output_name else input_path.stem
    output_file_path = (
        output_dir / f"{output_name}.{output_extension}"
    )  # Correct concatenation
    check_directory_exists(str(output_dir))  # Ensure output directory exists
    return str(output_file_path)


def add_identifier_name(
    name: str,
    default_identifier: Optional[str] = None,
    identifier: Optional[str] = None,
    suffix: Optional[str] = None,
) -> str:
    parts = [name.lower()]
    if default_identifier:
        parts.append(default_identifier)
    if identifier:
        parts.append(identifier)
    if suffix:
        parts.append(suffix)
    return "_".join(parts)


def add_suffix_to_filename(
    file_path: str, suffix: str, return_full_path: bool = True
) -> str:
    """
    Adds a suffix to the base name of a file, preserving the directory and extension.

    :param file_path: The original file path (e.g., "path/to/file.extension").
    :param suffix: The suffix to add to the base name (e.g., "_value").
    :param return_full_path: If True, returns the full path; if False, returns only the base name.
    :return: Modified file path or base name with the suffix.
    """
    path = Path(file_path)
    new_name = f"{path.stem}_{suffix}{path.suffix}"

    if return_full_path:
        return str(path.with_name(new_name))
    else:
        return new_name


def create_temp_directory(base_dir: str, sub_dir: str) -> str:
    temp_dir = os.path.join(base_dir, sub_dir)
    os.makedirs(temp_dir, exist_ok=True)
    return temp_dir


def cleanup_directory(directory: str) -> None:
    if os.path.exists(directory):
        shutil.rmtree(directory)
