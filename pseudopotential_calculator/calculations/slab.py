from __future__ import annotations

from typing import TYPE_CHECKING, Self, cast

import numpy as np
from ase import Atoms
from ase.build import (
    add_vacuum,  # type: ignore bad library
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
    calculator.param.spin_polarized = parameters.spin_polarized
    calculator.cell.kpoint_mp_grid = parameters.kpoint_mp_grid

    # Prevent the bulk cell from rotating
    calculator.cell.cell_constraints = "0 0 0\n0 0 0"
    calculator.task = "Energy"

    atoms.set_constraint(FixAtoms(mask=[True for _ in atoms]))  # type: ignore unknown
    calculator.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculator
    return calculator


def get_slab_relax_calculator(
    atoms: Atoms,
    n_fixed_layer: int,
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
    mask = [i < n_fixed_layer for i in range(len(atoms))]
    atoms.set_constraint(FixAtoms(mask))  # type: ignore unknown
    calculator.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculator
    return calculator


def _get_height_per_layer(
    atom: Atoms,
    slab_direction: tuple[int, int, int],
) -> int:
    # create a slab with two layers,
    # the hight per layer is the difference between the two heights
    slab_single_layer = cast(
        Atoms,
        surface(
            atom,
            slab_direction,
            layers=2,
            vacuum=0,
        ),
    )

    position_0 = slab_single_layer.positions[0, 2]  # type: ignore unknown
    position_1 = slab_single_layer.positions[1, 2]  # type: ignore unknown
    return np.abs(cast(float, position_0) - cast(float, position_1))


def get_surface(
    atom: Atoms,
    slab_direction: tuple[int, int, int],
    *,
    n_layer: int,
    n_vacuum_layer: int,
) -> Atoms:
    slab = cast(Atoms, surface(atom, slab_direction, n_layer, 0))

    height_per_layer = _get_height_per_layer(atom, slab_direction)
    add_vacuum(slab, n_vacuum_layer * height_per_layer)
    return slab


def _get_n_vacuum_layers_from_slab(atom: Atoms) -> float:
    cell_height = atom.get_cell().array[-1][-1]

    atom_heights = list[float](atom.positions[:, 2])  # type: ignore unknown
    slab_height = max(atom_heights) - min(atom_heights)
    height_per_layer = slab_height / len(atom_heights)

    vacuum_height = cell_height - slab_height

    return vacuum_height / height_per_layer


def plot_energy_against_n_vacuum_layer(
    calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    """Plot energy against number of vacuum layer.

    This assumes the vacuum is layered up in z direction.
    """
    energies = list[float]()
    n_vacuum_layer = list[float]()
    for calculator in calculators:
        atom = get_calculator_atom(calculator)

        energies.append(get_atom_potential_energy(atom))
        n_vacuum_layers = _get_n_vacuum_layers_from_slab(atom)
        n_vacuum_layer.append(n_vacuum_layers)

    return plot_data_comparison(
        ("Vacuum Height", np.array(n_vacuum_layer), "N layers"),
        ("Energy", np.array(energies), "J"),
        ax=ax,
    )


def get_atom_displacement(atom: Atoms, n_th_layer: int) -> float:
    heights = atom.positions[:, 2]  # type: ignore bad lib
    sorted_heights = np.unique(heights)[::-1]  # type: ignore bad lib
    return sorted_heights[n_th_layer - 1] - sorted_heights[n_th_layer]


def get_n_free_layer(atom: Atoms) -> int:
    heights = atom.positions[:, 2]  # type: ignore bad lib
    return (len(np.unique(heights)) + 1) / 2  # type: ignore bad lib


def plot_displacement_against_n_free_layer(
    calculators: list[Castep],
    n_th_layer: int,
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    """Plot displacement of free atoms against number of free layer."""
    displacement = list[float]()
    n_free_layer = list[int]()
    for calculator in calculators:
        atom = get_calculator_atom(calculator)

        displacement.append(
            get_atom_displacement(atom, n_th_layer),
        )
        n_free_layers = get_n_free_layer(atom)
        n_free_layer.append(n_free_layers)

    return plot_data_comparison(
        ("N Free Layers", np.array(n_free_layer), None),
        ("Displacement", np.array(displacement), "m"),
        ax=ax,
    )
