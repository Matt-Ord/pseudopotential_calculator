from __future__ import annotations

import os
import subprocess
import warnings
from dataclasses import dataclass, field
from functools import cache
from pathlib import PosixPath
from typing import TYPE_CHECKING

from pseudopotential_calculator.castep import (
    get_calculator_directory,
    get_calculator_label,
    prepare_calculator,
)

if TYPE_CHECKING:
    from pathlib import Path

    from ase.calculators.castep import Castep


@cache
def _get_hpc_account() -> str:
    try:
        return os.environ["HPC_ACCOUNT"]
    except KeyError:
        return input("HPC Account: ")


@dataclass
class HPCTaskConfig:
    """Configuration for a HPC task."""

    account: str = field(default_factory=_get_hpc_account)
    partition: str = "icelake-himem"
    n_nodes: int = field(default=1, kw_only=True)
    n_tasks: int = field(default=72, kw_only=True)
    time: str = field(default="12:00:00", kw_only=True)


def get_submit_script(calculator: Castep, config: HPCTaskConfig | None = None) -> str:
    config = HPCTaskConfig() if config is None else config
    label = get_calculator_label(calculator)
    return f"""#!/bin/bash
#SBATCH --account={config.account}
#SBATCH --partition={config.partition}
#SBATCH --job-name={label}
#SBATCH --nodes={config.n_nodes}
#SBATCH --ntasks={config.n_tasks}
#SBATCH --time={config.time}

. /etc/profile.d/modules.sh


module purge
module load rhel8/default-icl castep/impi/22.11

mpirun castep.mpi {label}
"""


def _get_submit_file_path(calculator: Castep) -> Path:
    directory = get_calculator_directory(calculator)
    return directory / "submit.sh"


def prepare_submit_script_for_calculator(
    calculator: Castep,
    config: HPCTaskConfig | None = None,
) -> None:
    file_path = _get_submit_file_path(calculator)

    script = get_submit_script(calculator, config)
    file_path.write_text(script)


def prepare_calculator_with_submit_script(
    calculator: Castep,
    config: HPCTaskConfig | None = None,
) -> None:
    prepare_calculator(calculator)
    prepare_submit_script_for_calculator(calculator, config)


def prepare_all_submit_scripts(
    calculations: list[Castep],
    config: HPCTaskConfig,
) -> None:
    for calculation in calculations:
        prepare_submit_script_for_calculator(calculation, config)


def _get_submit_all_script_for_directories(
    directories: list[Path],
) -> str:
    directories_string = " ".join(f'"{d}"' for d in directories)
    return f"""#!/bin/bash

#!/bin/bash

directories=({directories_string})

for i in "${{directories[@]}}"; do
    # Navigate to the directory
    cd "${{i}}" || {{ echo "Failed to cd into directory"; exit 1; }}

    # Submit the SLURM job
    sbatch submit.sh

    # Return to the previous directory
    cd - || {{ echo "Failed to cd back"; exit 1; }}

    # Print the status message
    echo "${{i}} done"
done
"""


def _get_submit_all_script(
    calculations: list[Castep],
    directory: Path,
) -> str:
    directories = [
        get_calculator_directory(calculation) for calculation in calculations
    ]
    directories = [d.relative_to(directory) for d in directories]
    return _get_submit_all_script_for_directories(directories)


def _try_grant_execute_permissions_to_file(path: Path) -> None:
    """Grant execute permissions to all files in the given directory."""
    # Construct the chmod command
    command = ["chmod", "-R", "+x", f"{path.absolute()}"]
    try:
        subprocess.run(command, check=True, capture_output=True)  # noqa: S603
    except subprocess.CalledProcessError:
        warnings.warn(
            f"Unable to make {path.absolute()} excecutable, "
            f"to do this manually call `chmod -R +x {path.absolute()}`",
            stacklevel=1,
        )


def prepare_submit_all_script(calculations: list[Castep], directory: Path) -> None:
    script = _get_submit_all_script(calculations, directory)

    file_path = directory / "submit.sh"
    file_path.write_text(script)
    _try_grant_execute_permissions_to_file(file_path)


def _get_ssh_username() -> str:
    try:
        return os.environ["HPC_USERNAME"]
    except KeyError:
        return input("Username: ")


def _get_hpc_workspace_directory() -> PosixPath:
    try:
        return PosixPath(os.environ["HPC_WORKSPACE"])
    except KeyError:
        username = _get_ssh_username()
        return PosixPath(f"/rds/user/{username}/hpc-work")


def _get_relative_hpc_path(remote_folder: PosixPath) -> PosixPath:
    workspace_folder = _get_hpc_workspace_directory()
    return workspace_folder / remote_folder


def copy_files_to_hpc(local_folder: Path, remote_folder: PosixPath) -> None:
    """Make use of the scp command to copy files from local to remote."""
    username = _get_ssh_username()

    remote_folder_absolute = _get_relative_hpc_path(remote_folder)
    command = [
        "scp",
        "-r",
        f"{local_folder.absolute()}",
        f"{username}@login.hpc.cam.ac.uk:{remote_folder_absolute.as_posix()}",
    ]
    subprocess.run(
        command,  # noqa: S603
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        check=True,
    )


def copy_files_from_hpc(local_folder: Path, remote_folder: PosixPath) -> None:
    """Make use of the scp command to copy files from local to remote."""
    username = _get_ssh_username()

    remote_folder_absolute = _get_relative_hpc_path(remote_folder)
    command = [
        "scp",
        "-r",
        f"{username}@login.hpc.cam.ac.uk:{remote_folder_absolute.as_posix()}",
        f"{local_folder.absolute()}",
    ]
    subprocess.run(
        command,  # noqa: S603
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        check=True,
    )
