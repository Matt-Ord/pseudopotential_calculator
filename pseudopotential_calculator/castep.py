from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Self, cast

from ase import Atoms
from ase.calculators.castep import (
    Castep,
    _get_indices_to_sort_back,  # type: ignore private
)
from ase.io.castep import (
    read_bands,  # type: ignore unknown
    read_castep_castep,  # type: ignore unknown
    read_castep_cell,  # type: ignore unknown
)


@dataclass
class CastepConfig:
    """Configuration for a CASTEP calculation."""

    directory: Path
    label: str


def get_calculator_directory(calculator: Castep) -> Path:
    return Path(calculator._directory)  # type: ignore this is the only way # noqa: SLF001


def get_calculator_atom(calculator: Castep) -> Atoms:
    return calculator.atoms  # type: ignore this is the only way


def get_calculator_label(calculator: Castep) -> str:
    return calculator._label  # type: ignore this is the only way # noqa: SLF001


def get_calculator_config(calculator: Castep) -> CastepConfig:
    return CastepConfig(
        get_calculator_directory(calculator),
        get_calculator_label(calculator),
    )


# Fixes bug reported here <https://gitlab.com/ase/ase/-/issues/1531>
class _CastepFixed(Castep):
    def read(self: Self, castep_file: str) -> None:
        """Read a castep file into the current instance."""
        atoms = cast(Atoms, read_castep_castep(castep_file))

        self.results = atoms.calc.results  # type: ignore a

        self._cut_off_energy = atoms.calc._cut_off_energy  # type: ignore a  # noqa: SLF001
        for k, v in atoms.calc._parameters_header.items():  # type: ignore a  # noqa: SLF001
            setattr(self.param, k, v)  # type: ignore a

        if self.atoms and not self._set_atoms:  # type: ignore a
            # compensate for internal reordering of atoms by CASTEP
            # using the fact that the order is kept within each species

            indices = _get_indices_to_sort_back(
                self.atoms.symbols,  # type: ignore a
                atoms.symbols,
            )
            positions_frac_atoms = atoms.get_scaled_positions()[indices]  # type: ignore a
            self.atoms.set_scaled_positions(positions_frac_atoms)  # type: ignore a
            keys = [
                "forces",
                "charges",
                "magmoms",
                "hirshfeld_volume_ratios",
                "hirshfeld_charges",
                "hirshfeld_magmoms",
            ]
            for k in keys:
                if k not in self.results:  # type: ignore a
                    continue
                self.results[k] = self.results[k][indices]  # type: ignore a

        else:
            # ------------ modified from upstream -------------- #
            self._kpoints = atoms.calc._kpoints  # type: ignore a  # noqa: SLF001
            self._species_pot = atoms.calc._species_pot  # type: ignore a  # noqa: SLF001
            self.cell.species_pot = atoms.calc._species_pot  # type: ignore a  # noqa: SLF001
            self._total_time = atoms.calc._total_time  # type: ignore a  # noqa: SLF001
            self._peak_memory = atoms.calc._peak_memory  # type: ignore a  # noqa: SLF001
            # -------------------------------------------------- #

            atoms.set_initial_charges(self.results.get("charges"))  # type: ignore a
            atoms.set_initial_magnetic_moments(self.results.get("magmoms"))  # type: ignore a
            atoms.calc = self

        self._kpoints = atoms.calc._kpoints  # type: ignore a  # noqa: SLF001

        self.cell.species_pot = atoms.calc._species_pot  # type: ignore a  # noqa: SLF001

        self._total_time = atoms.calc._total_time  # type: ignore a  # noqa: SLF001
        self._peak_memory = atoms.calc._peak_memory  # type: ignore a  # noqa: SLF001

        # Read in eigenvalues from bands file
        bands_file = castep_file[:-7] + ".bands"
        if (
            self.param.task.value is not None
            and self.param.task.value.lower() == "bandstructure"
        ):
            self._band_structure = self.band_structure(bandfile=bands_file)  # type: ignore a
        else:
            try:
                (
                    self._ibz_kpts,
                    self._ibz_weights,
                    self._eigenvalues,  # type: ignore a
                    self._efermi,  # type: ignore a
                ) = read_bands(bands_file)
            except FileNotFoundError:
                warnings.warn(  # noqa: B028
                    "Could not load .bands file, eigenvalues and "
                    "Fermi energy are unknown",
                )


def get_default_calculator(
    config: CastepConfig,
) -> Castep:
    """Set up the CASTEP calculation parameters using the values from config."""
    calculator = _CastepFixed(
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


def load_calculator(config: CastepConfig) -> Castep:
    calculator = get_default_calculator(config)

    castep_file = f"{config.label}.castep"
    castep_path = config.directory / castep_file
    calculator.read(f"{castep_path}")  # type: ignore unknown

    cell_file = f"{config.label}.cell"
    cell_path = config.directory / cell_file
    with cell_path.open() as f:
        cell_atoms = cast(Atoms, read_castep_cell(f))
        calculator.cell = cell_atoms.calc.cell  # type: ignore bad lib type

    calculator.push_oldstate()
    return calculator


def load_all_calculators(directory: Path) -> list[Castep]:
    out = list[Castep]()
    for root, _, files in os.walk(directory):
        if any(f.endswith(".err") for f in files):
            warnings.warn(f".err file in {root}, skipping", stacklevel=1)
            continue

        for file in files:
            if file.endswith(".castep"):
                root_directory = Path(root)
                config = CastepConfig(root_directory, file.removesuffix(".castep"))
                out.append(load_calculator(config))

    return out
