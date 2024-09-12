from __future__ import annotations

from typing import TYPE_CHECKING, Literal, cast, overload

import numpy as np

from pseudopotential_calculator.castep import (
    get_atom_potential_energy,
    get_calculator_atom,
    get_calculator_cutoff_energy,
)
from pseudopotential_calculator.util import plot_data_comparison

if TYPE_CHECKING:
    from ase.calculators.castep import Castep
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

XCFunctional = Literal["PBE", "LDA", "WC"]
Task = Literal["Energy", "GeometryOptimization"]


def _get_cell_lengths_from_calculator(
    calculator: Castep,
) -> tuple[float, ...]:
    atom = get_calculator_atom(calculator)
    return tuple(atom.get_cell_lengths_and_angles()[0:3])


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


def plot_theoretical_cell_length(ax: Axes, length: float) -> None:
    ax.axhline(y=length, linestyle="--", label=f"Theoretical Length: {length} Å")  # type: ignore bad library


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

    return plot_data_comparison(
        ("N K Points", np.array(n_k_points), None),
        ("Bond Length", np.array(bond_lengths), "m"),
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
        energies.append(get_calculator_atom(calculator).get_potential_energy())  # type: ignore inkown
        n_k_points.append(_get_n_k_points_from_calculator(calculator, direction))

    return plot_data_comparison(
        ("N K Points", np.array(n_k_points), None),
        ("Energy", np.array(energies), "J"),
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
        atom = get_calculator_atom(calculator)

        energies.append(get_atom_potential_energy(atom))
        cutoff_energy.append(get_calculator_cutoff_energy(calculator))

    return plot_data_comparison(
        ("Cutoff Energy", np.array(cutoff_energy), "J"),
        ("Energy", np.array(energies), "J"),
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
        cutoff_energy.append(get_calculator_cutoff_energy(calculator))

    return plot_data_comparison(
        ("Cutoff Energy", np.array(cutoff_energy), "J"),
        ("Cell Length", np.array(cell_length), r"$/AA$"),
        ax=ax,
    )
