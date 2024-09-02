from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, NamedTuple, Self, cast, overload

import numpy as np
from ase import Atoms

from pseudopotential_calculator.calculations.generic import OptimizationParamsBase
from pseudopotential_calculator.castep import (
    CastepConfig,
    get_calculator_atom,
    get_default_calculator,
)
from pseudopotential_calculator.util import plot_data_comparison

if TYPE_CHECKING:
    from ase.calculators.castep import Castep
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

XCFunctional = Literal["PBE", "LDA", "WC"]


class BulkOptimizationParams(OptimizationParamsBase):
    """Parameters of a bulk calculation."""

    @property
    def kpoint_mp_grid(self: Self) -> str:
        return f"{self.n_k_points} {self.n_k_points} {self.n_k_points}"


def get_bulk_optimization_calculator(
    atoms: Atoms,
    parameters: BulkOptimizationParams,
    config: CastepConfig,
) -> Castep:
    calculator = get_default_calculator(config)

    calculator.param.xc_functional = parameters.xc_functional
    calculator.param.cut_off_energy = parameters.cut_off_energy
    calculator.param.spin_polarized = parameters.spin_polarized
    calculator.cell.kpoint_mp_grid = parameters.kpoint_mp_grid

    # Prevent the bulk cell from rotating
    calculator.cell.cell_constraints = "1 1 1\n0 0 0"
    calculator.task = "GeometryOptimization"

    calculator.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculator
    return calculator


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


def _get_cutoff_energy_from_calculator(
    calculator: Castep,
) -> float:
    return cast(
        float,
        calculator.param.cut_off_energy.raw_value[0],  # type: ignore unknown
    )


def plot_theoretical_cell_length(ax: Axes) -> None:
    ax.axhline(y=2.53, linestyle="--", label="Theoretical Length: 2.53")  # type: ignore bad library


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
        n_k_points.append(_get_n_k_points_from_calculator(calculator, direction))

    class PlotTuple(NamedTuple):
        n_k_points: np.ndarray[Any, np.dtype[np.float64]]
        bond_lengths: np.ndarray[Any, np.dtype[np.float64]]

    p = PlotTuple(np.array(n_k_points), np.array(bond_lengths))
    return plot_data_comparison(
        p,
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

    class PlotTuple(NamedTuple):
        n_k_points: np.ndarray[Any, np.dtype[np.float64]]
        energies: np.ndarray[Any, np.dtype[np.float64]]

    p = PlotTuple(np.array(n_k_points), np.array(energies))
    return plot_data_comparison(
        p,
        ax=ax,
    )


def plot_energy_against_cutoff_energy(
    calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    cutoff_energy = list[float]()
    energies = list[float]()
    for calculator in calculators:
        energies.append(
            cast(Atoms, get_calculator_atom(calculator)).get_potential_energy(),  # type: ignore
        )  # type: ignore inkown
        cutoff_energy.append(_get_cutoff_energy_from_calculator(calculator))

    class PlotTuple(NamedTuple):
        cutoff_energy: np.ndarray[Any, np.dtype[np.float64]]
        energies: np.ndarray[Any, np.dtype[np.float64]]

    p = PlotTuple(np.array(cutoff_energy), np.array(energies))
    return plot_data_comparison(
        p,
        ax=ax,
    )


def plot_cell_length_against_cutoff_energy(
    calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    cutoff_energy = list[float]()
    cell_length = list[float]()
    for calculator in calculators:
        cell_length.append(_get_cell_lengths_from_calculator(calculator)[0])
        cutoff_energy.append(_get_cutoff_energy_from_calculator(calculator))

    class PlotTuple(NamedTuple):
        cutoff_energy: np.ndarray[Any, np.dtype[np.float64]]
        cell_length: np.ndarray[Any, np.dtype[np.float64]]

    p = PlotTuple(np.array(cutoff_energy), np.array(cell_length))
    return plot_data_comparison(
        p,
        ax=ax,
    )
