from __future__ import annotations

import os
from dataclasses import dataclass, field
from functools import cache
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
