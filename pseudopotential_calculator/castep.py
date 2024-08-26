from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from ase.calculators.castep import Castep


@dataclass
class CastepConfig:
    """Configuration for a CASTEP calculation."""

    directory: Path
    label: str


def get_calculator_directory(calculator: Castep) -> Path:
    return Path(calculator._directory)  # type: ignore this is the only way # noqa: SLF001


def get_calculator_label(calculator: Castep) -> str:
    return calculator._label  # type: ignore this is the only way # noqa: SLF001


def get_calculator_config(calculator: Castep) -> CastepConfig:
    return CastepConfig(
        get_calculator_directory(calculator),
        get_calculator_label(calculator),
    )


def get_default_calculator(
    config: CastepConfig,
) -> Castep:
    """Set up the CASTEP calculation parameters using the values from config."""
    calculator = Castep(
        directory=config.directory.as_posix(),
        label=config.label,
        keyword_tolerance=0,
    )

    # Don't prefix calculations with 0000
    calculator._track_output = False  # type: ignore only way # noqa: SLF001
    calculator._try_reuse = True  # type: ignore only way # noqa: SLF001
    calculator._pedantic = True  # type: ignore only way # noqa: SLF001
    calculator._set_atoms = True  # type: ignore only way # noqa: SLF001
    calculator.param.num_dump_cycles = 0
    calculator.param.reuse = True

    return calculator


def prepare_calculator(calculator: Castep) -> None:
    calculator.prepare_input_files()  # type: ignore bad lib types


def load_all_calculators(directory: Path) -> list[Castep]:
    out = list[Castep]()
    for root, _, files in os.walk(directory):
        for file in files:
            if Path(file).suffix == ".castep":
                root_path = Path(root)
                path = f"{(Path(root) / file).absolute()}"
                calculator = get_default_calculator(CastepConfig(root_path, file))
                calculator.read(path)  # type: ignore unknown
                out.append(calculator)

    return out
