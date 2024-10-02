from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Self, cast

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


@dataclass
class CastepParams:
    """class for castep parameters."""

    symmetry_generate: bool = field(default=False, kw_only=True)
    n_k_points: tuple[int, int, int] = field(default=(1, 1, 1), kw_only=True)
    cut_off_energy: float = field(default=600, kw_only=True)
    xc_functional: str = field(
        default="PBE",
        kw_only=True,
    )
    spin_polarized: bool = field(default=False, kw_only=True)
    max_scf_cycles: int = 30
    length_constraints: tuple[int, int, int] = field(default=(0, 0, 0), kw_only=True)
    angle_constraints: tuple[int, int, int] = field(default=(0, 0, 0), kw_only=True)
    continuation: str = field(
        default="default",
        kw_only=True,
    )

    @property
    def cell_constraints(self: Self) -> str:
        return (
            f"{self.length_constraints[0]} {self.length_constraints[1]} "
            f"{self.length_constraints[2]}\n"
            f"{self.angle_constraints[0]} {self.angle_constraints[1]} "
            f"{self.angle_constraints[2]}"
        )


def get_calculator(
    atoms: Atoms,
    parameters: CastepParams,
    config: CastepConfig,
    task: Literal["Energy", "GeometryOptimization"] = "GeometryOptimization",
) -> Castep:
    calculator = get_default_calculator(config)
    calculator.task = task

    calculator.param.xc_functional = parameters.xc_functional
    calculator.param.cut_off_energy = parameters.cut_off_energy
    calculator.param.spin_polarized = parameters.spin_polarized
    calculator.param.max_scf_cycles = parameters.max_scf_cycles
    calculator.param.continuation = parameters.continuation
    calculator.cell.kpoint_mp_grid = (
        f"{parameters.n_k_points[0]}"
        f"{parameters.n_k_points[1]} "
        f"{parameters.n_k_points[2]}"
    )
    calculator.cell.symmetry_generate = parameters.symmetry_generate
    calculator.cell.cell_constraints = parameters.cell_constraints

    calculator.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculator
    return calculator


def get_calculator_directory(calculator: Castep) -> Path:
    return Path(calculator._directory)  # type: ignore this is the only way # noqa: SLF001


def get_calculator_atom(calculator: Castep) -> Atoms:
    atoms = cast(Atoms | None, calculator.atoms)  # type: ignore this is the only way
    if atoms is None:
        msg = "Calculator did not have Atoms attached."
        raise ValueError(msg)
    return atoms


def get_calculator_label(calculator: Castep) -> str:
    return calculator._label  # type: ignore this is the only way # noqa: SLF001


def get_calculator_config(calculator: Castep) -> CastepConfig:
    return CastepConfig(
        get_calculator_directory(calculator),
        get_calculator_label(calculator),
    )


def get_atom_potential_energy(atom: Atoms) -> float:
    return atom.get_potential_energy()  # type: ignore this is the only way


def get_calculator_cutoff_energy(
    calculator: Castep,
) -> float:
    return cast(float, calculator.param.cut_off_energy.raw_value[0])  # type: ignore unknown


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
    calculator._pedantic = True  # type: ignore only way # noqa: SLF001
    calculator._set_atoms = True  # type: ignore only way # noqa: SLF001
    calculator.param.num_dump_cycles = 0
    calculator.param.elec_energy_tol = 1e-06
    calculator.param.geom_energy_tol = 1e-05
    calculator.param.geom_disp_tol = 1e-03

    return calculator


def prepare_calculator(calculator: Castep) -> None:
    calculator.prepare_input_files()  # type: ignore bad lib types


def get_calculator_old_atoms(calculator: Castep) -> Atoms:
    return cast(Atoms, calculator.__old_atoms)  # type: ignore bad lib # noqa: SLF001


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
        calculator.__old_atoms = cell_atoms  # type: ignore modified originally not useful _old_atoms# noqa: SLF001

    calculator.push_oldstate()
    return calculator


def load_calculator_atoms(config: CastepConfig) -> Atoms:
    return get_calculator_atom(load_calculator(config))


def load_all_calculators(directory: Path) -> list[Castep]:
    out = list[Castep]()
    for root, _, files in os.walk(directory):
        if any(f.endswith(".err") for f in files):
            warnings.warn(f".err file in {root}, skipping", stacklevel=1)
            continue

        for file in files:
            if file.endswith(".castep"):
                try:
                    root_directory = Path(root)
                    config = CastepConfig(root_directory, file.removesuffix(".castep"))
                    out.append(load_calculator(config))
                except AttributeError:
                    warnings.warn(
                        f"unable to load calculator in {root}, skipping",
                        stacklevel=1,
                    )

    return out
