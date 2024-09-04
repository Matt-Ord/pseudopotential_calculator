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
    get_atom_potential_energy,
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


def get_slab_vacuum_calculator(
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


def get_heigh_per_vacuum_layer(
    atom: Atoms,
    slab_direction: tuple[int, int, int],
) -> int:
    slab_single_layer = cast(
        Atoms,
        surface(
            atom,
            slab_direction,
            layers=2,
            vacuum=0,
        ),
    )

    return slab_single_layer.get_cell().array[-1][-1]


def get_surface(
    atom: Atoms,
    slab_direction: tuple[int, int, int],
    n_layer: int,
    n_vacuum_layer: int,
) -> Atoms:
    height_per_vacuum_layer = get_heigh_per_vacuum_layer(atom, slab_direction)
    slab = cast(
        Atoms,
        surface(
            atom,
            slab_direction,
            n_layer,
            n_vacuum_layer * height_per_vacuum_layer / 2,
        ),
    )  # type: ignore bad lib
    zmin = cast(float, np.min(slab.positions[:, 2]))  # type: ignore bad lib
    # Adjust the positions to move the bottom layer to z = 0
    slab.translate([0, 0, -zmin])  # type: ignore bad lib
    return slab


class _PlotTuple(NamedTuple):
    n_vacuum_layer: np.ndarray[Any, np.dtype[np.float64]]
    energies: np.ndarray[Any, np.dtype[np.float64]]


def get_n_vacuum_layers_from_out_put(atom: Atoms) -> int:
    cell_height = atom.get_cell().array[-1][-1]
    slab_height = atom.get_cell().array[-1][-1]
    slab_height_minus_1_layer = cast(float, atom.get_positions()[-2][-1])  # type: ignore bad lib
    height_per_vacuum_layer = slab_height - slab_height_minus_1_layer
    return round((cell_height - slab_height) / height_per_vacuum_layer)


def plot_energy_against_n_vacuum_layer(
    calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    """Plot energy against number of vacuum layer.

    This assumes the vacuum is layered up in z direction.
    """
    energies = list[float]()
    n_vacuum_layer = list[int]()
    for calculator in calculators:
        atom = get_calculator_atom(calculator)
        energies.append(
            get_atom_potential_energy(atom),
        )
        n_vacuum_layers = get_n_vacuum_layers_from_out_put(atom)
        n_vacuum_layer.append(n_vacuum_layers)

    p = _PlotTuple(np.array(n_vacuum_layer), np.array(energies))
    return plot_data_comparison(
        p,
        ax=ax,
    )
