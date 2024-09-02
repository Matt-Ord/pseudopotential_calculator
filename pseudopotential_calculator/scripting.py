import os
from pathlib import Path, PosixPath
from typing import Literal

from pseudopotential_calculator.hpc import copy_files_from_hpc, copy_files_to_hpc


def prompt_for_decision(question: str, default: Literal["y", "n"] = "y") -> bool:
    valid_responses = {"y": True, "n": False, "yes": True, "no": False}

    if default not in ["y", "n"]:
        msg = f"Invalid default answer: '{default}'"
        raise TypeError(msg)

    prompt_msg = {"y": " [Y/n]: ", "n": " [y/N]: "}[default]

    while True:
        response = input(f"{question} {prompt_msg}").strip().lower()

        if response == "":
            return valid_responses[default]
        if response in valid_responses:
            return valid_responses[response]
        print("Please respond with 'y' or 'n'.")  # noqa: T201


def maybe_copy_files_from_hpc(local_folder: Path, remote_folder: PosixPath) -> None:
    """Make use of the scp command to copy files from local to remote."""
    skip = os.environ.get("SKIP_DOWNLOAD", False)
    if skip or not prompt_for_decision("Copy files from HPC"):
        return
    copy_files_from_hpc(local_folder, remote_folder)


def maybe_copy_files_to_hpc(local_folder: Path, remote_folder: PosixPath) -> None:
    """Make use of the scp command to copy files to local to remote."""
    skip = os.environ.get("SKIP_UPLOAD", False)
    if skip or not prompt_for_decision("Copy files to HPC"):
        return
    copy_files_to_hpc(local_folder, remote_folder)
