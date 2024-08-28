from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Self, cast, overload

import numpy as np
from ase import Atoms

from pseudopotential_calculator.castep import CastepConfig, get_default_calculator
from pseudopotential_calculator.util import get_figure

if TYPE_CHECKING:
    from ase.calculators.castep import Castep
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D


@dataclass
class BulkOptimizationParams:
    """Parameters of a bulk calculation."""

    n_k_points: int = field(default=1, kw_only=True)
    cut_off_energy: float = field(default=340, kw_only=True)
    xc_functional: str = field(default="PBE", kw_only=True)

    @property
    def kpoint_mp_grid(self: Self) -> str:
        return f"{self.n_k_points} {self.n_k_points} {self.n_k_points}"


def get_bulk_optimization_calculator(
    atoms: Atoms,
    parameters: BulkOptimizationParams,
    config: CastepConfig,
) -> Castep:
    calculation = get_default_calculator(config)

    calculation.param.xc_functional = parameters.xc_functional
    calculation.param.cut_off_energy = parameters.cut_off_energy

    calculation.cell.kpoint_mp_grid = parameters.kpoint_mp_grid
    # Prevent the bulk cell from rotating
    calculation.cell.cell_constraints = "1 1 1\n0 0 0"
    calculation.task = "GeometryOptimization"

    calculation.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculation
    return calculation


def _plot_cell_length_against_n_k_points(
    n_k_points: np.ndarray[Any, np.dtype[np.float64]],
    bond_lengths: np.ndarray[Any, np.dtype[np.float64]],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    fig, ax = get_figure(ax)

    args = np.argsort(n_k_points)
    (line,) = ax.plot(n_k_points[args], bond_lengths[args])  # type: ignore library
    ax.set_xlabel("number of k-points")  # type: ignore library
    ax.set_ylabel(r"Cell Length / $\AA$")  # type: ignore library
    ax.set_title("Plot of cell length vs number of k-points")  # type: ignore library
    return fig, ax, line


def _plot_energy_against_n_k_points(
    n_k_points: np.ndarray[Any, np.dtype[np.float64]],
    energy: np.ndarray[Any, np.dtype[np.float64]],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    fig, ax = get_figure(ax)

    args = np.argsort(n_k_points)
    (line,) = ax.plot(n_k_points[args], energy[args])  # type: ignore library
    ax.set_xlabel("number of k-points")  # type: ignore library
    ax.set_ylabel("Energy / eV")  # type: ignore library
    ax.set_title("Plot of energy vs number of k-points")  # type: ignore library
    return fig, ax, line


def _get_cell_lengths_from_calculator(
    calculator: Castep,
) -> tuple[float, ...]:
    return cast(Atoms, calculator.atoms).get_cell_lengths_and_angles()[0:3]  # type: ignore unknown


@overload
def _get_n_k_points_from_calculator(
    calculator: Castep,
    direction: None = None,
) -> tuple[int, ...]:
    ...


@overload
def _get_n_k_points_from_calculator(
    calculator: Castep,
    direction: int,
) -> int:
    ...


def _get_n_k_points_from_calculator(
    calculator: Castep,
    direction: int | None = None,
) -> tuple[int, ...] | int:
    kpoint_mp_grid = cast(
        tuple[int, ...],
        calculator.cell.kpoint_mp_grid.raw_value,  # type: ignore unknown
    )
    if direction is None:
        return kpoint_mp_grid
    return kpoint_mp_grid[direction]


def plot_cell_length_against_n_k_points(
    calculators: list[Castep],
    direction: int = 0,
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    n_k_points = list[float]()
    bond_lengths = list[float]()
    for calculator in calculators:
        bond_lengths.append(_get_cell_lengths_from_calculator(calculator)[direction])
        print(_get_cell_lengths_from_calculator(calculator)[direction])
        n_k_points.append(_get_n_k_points_from_calculator(calculator, direction))

    return _plot_cell_length_against_n_k_points(
        np.array(n_k_points),
        np.array(bond_lengths),
        ax=ax,
    )


def plot_energy_against_n_k_points(
    calculators: list[Castep],
    direction: int = 0,
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    n_k_points = list[float]()
    energies = list[float]()
    for calculator in calculators:
        energies.append(cast(Atoms, calculator.atoms).get_potential_energy())  # type: ignore inkown
        n_k_points.append(_get_n_k_points_from_calculator(calculator, direction))

    return _plot_energy_against_n_k_points(
        np.array(n_k_points),
        np.array(energies),
        ax=ax,
    )


def _get_cutoff_energy_from_calculator(
    calculator: Castep,
) -> float:
    return cast(
        float,
        calculator.param.cut_off_energy.raw_value[0],  # type: ignore unknown
    )


def _plot_energy_against_cutoff_energy(
    cutoff_energy: np.ndarray[Any, np.dtype[np.float64]],
    energy: np.ndarray[Any, np.dtype[np.float64]],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    fig, ax = get_figure(ax)

    args = np.argsort(cutoff_energy)
    (line,) = ax.plot(cutoff_energy[args], energy[args])  # type: ignore library
    ax.set_xlabel("Cuttoff energy /eV")  # type: ignore library
    ax.set_ylabel("Energy / eV")  # type: ignore library
    ax.set_title("Plot of energy vs cutoff energy")  # type: ignore library
    return fig, ax, line


def plot_energy_against_cutoff_energy(
    calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    cutoff_energy = list[float]()
    energies = list[float]()
    for calculator in calculators:
        energies.append(cast(Atoms, calculator.atoms).get_potential_energy())  # type: ignore inkown
        print(cast(Atoms, calculator.atoms).get_potential_energy())
        cutoff_energy.append(_get_cutoff_energy_from_calculator(calculator))

    return _plot_energy_against_cutoff_energy(
        np.array(cutoff_energy),
        np.array(energies),
        ax=ax,
    )
