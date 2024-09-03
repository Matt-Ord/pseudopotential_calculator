from __future__ import annotations

from typing import TYPE_CHECKING, Any, NamedTuple, Self, cast

import numpy as np
from ase import Atoms
from ase.build import (
    surface,  # type: ignore bad library
)
from ase.constraints import FixAtoms

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


class SlabOptimizationParams(OptimizationParamsBase):
    """Parameters of a slab calculation."""

    @property
    def kpoint_mp_grid(self: Self) -> str:
        return f"{self.n_k_points} {self.n_k_points} {1}"


def get_slab_vaccum_calculator(
    atoms: Atoms,
    parameters: SlabOptimizationParams,
    config: CastepConfig,
) -> Castep:
    calculator = get_default_calculator(config)

    calculator.param.xc_functional = parameters.xc_functional
    calculator.param.cut_off_energy = parameters.cut_off_energy
    calculator.param.spin_polarized = "true" if parameters.spin_polarized else "false"
    calculator.cell.kpoint_mp_grid = parameters.kpoint_mp_grid

    # Prevent the bulk cell from rotating
    calculator.cell.cell_constraints = "0 0 0\n0 0 0"
    calculator.task = "Energy"

    atoms.set_constraint(FixAtoms(mask=[True for _ in atoms]))  # type: ignore unknown
    calculator.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculator
    return calculator


def get_surface(
    atom: Atoms,
    slab_direction: tuple[int, int, int],
    n_layer: int,
    n_vaccum_layer: int,
) -> Atoms:
    # vaccum = distance between adjacent fcc 111 planes, calculated using geometry
    height_per_vaccum_layer = atom.cell.lengths()[0] * np.sqrt(2) / np.sqrt(3)

    slab = cast(
        Atoms,
        surface(
            atom,
            slab_direction,
            n_layer,
            n_vaccum_layer * height_per_vaccum_layer / 2,
        ),
    )  # type: ignore bad library
    zmin = np.min(slab.positions[:, 2])  # type: ignore bad library
    # Adjust the positions to move the bottom layer to z = 0
    slab.positions[:, 2] -= zmin  # type: ignore bad library

    return slab


def plot_energy_against_n_vaccum_layer(
    calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    """Plot energy against number of vaccum layer.

    This assumes the vaccum is layered up in z direction.
    """
    energies = list[float]()
    n_vaccum_layer = list[int]()
    for calculator in calculators:
        atom = get_calculator_atom(calculator)
        if atom is None:
            continue
        atom = cast(Atoms, get_calculator_atom(calculator))
        energies.append(
            atom.get_potential_energy(),  # type: ignore bad library
        )

        cell_height = atom.get_cell().array[-1][-1]
        slab_height = cast(float, atom.get_positions()[-1][-1])  # type: ignore
        slab_height_minus_1_layer = cast(float, atom.get_positions()[-2][-1])  # type: ignore

        height_per_vaccum_layer = slab_height - slab_height_minus_1_layer
        print(height_per_vaccum_layer)
        n_vaccum_layers = (cell_height - slab_height) / height_per_vaccum_layer
        n_vaccum_layer.append(n_vaccum_layers)

    class PlotTuple(NamedTuple):
        n_vaccum_layer: np.ndarray[Any, np.dtype[np.float64]]
        energies: np.ndarray[Any, np.dtype[np.float64]]

    p = PlotTuple(np.array(n_vaccum_layer), np.array(energies))
    return plot_data_comparison(
        p,
        ax=ax,
    )
