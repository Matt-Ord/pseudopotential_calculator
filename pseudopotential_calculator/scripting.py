import os
from pathlib import Path, PosixPath
from typing import Literal

from pseudopotential_calculator.hpc import copy_files_from_hpc


def prompt_for_decision(question: str, default: Literal["y", "n"] = "y") -> bool:
    valid_responses = {"y": True, "n": False, "yes": True, "no": False}

    if default not in valid_responses:
        msg = f"Invalid default answer: '{default}'"
        raise ValueError(msg)

    prompt_map = {
        "y": " [Y/n]: ",
        "n": " [y/N]: ",
    }

    while True:
        response = input(f"{question} {prompt_map[default]}").strip().lower()

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
