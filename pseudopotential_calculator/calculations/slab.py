from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, Self, cast

import numpy as np
from ase import Atoms
from ase.build import (
    add_vacuum,  # type: ignore bad library
    surface,  # type: ignore bad library
)
from ase.constraints import FixAtoms, FixedLine

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


def _get_heights_sorted(atom: Atoms) -> np.ndarray[Any, np.dtype[np.float64]]:
    heights = cast(np.ndarray[Any, np.dtype[np.float64]], atom.positions[:, 2])  # type: ignore bad lib
    return np.unique(heights)[::-1]


def get_n_th_slab_layer_height(atom: Atoms, n_th_layer: int) -> float:
    heights = _get_heights_sorted(atom)
    return heights[n_th_layer]


def get_height_per_slab_layer(atom: Atoms) -> float:
    _get_heights_sorted(atom)
    heights = _get_heights_sorted(atom)
    sorted_heights = np.unique(heights)[::-1]
    return sorted_heights[-2]


def get_n_slab_layer(atoms: Atoms) -> int:
    heights = _get_heights_sorted(atoms)
    return len(heights)


def get_slab_width(atoms: Atoms) -> int:
    # n_slab_atoms / n_atoms_per_repeat gives the number of slab units set to be repeated.
    # take square root to get width since the repeat s set to be the same in x and y directions
    n_slab_atoms = len(atoms.positions) - 1  # type: ignore bad lib
    n_atoms_per_repeat = get_n_slab_layer(atoms) - 1
    return np.sqrt(round(n_slab_atoms / n_atoms_per_repeat))


def plot_energy_against_slab_width(
    calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    width = list[int]()
    energies = list[float]()
    for calculator in calculators:
        atom = get_calculator_atom(calculator)
        energies.append(get_atom_potential_energy(atom))
        width.append(get_slab_width(atom))

    return plot_data_comparison(
        ("width", np.array(width), "atom_wide"),
        ("Energy", np.array(energies), "J"),
        ax=ax,
    )


def get_slab_calculator(
    atoms: Atoms,
    parameters: SlabOptimizationParams,
    config: CastepConfig,
    *,
    task: Literal["Energy", "GeometryOptimization"] = "Energy",
) -> Castep:
    calculator = get_default_calculator(config)

    calculator.param.xc_functional = parameters.xc_functional
    calculator.param.cut_off_energy = parameters.cut_off_energy
    calculator.param.spin_polarized = parameters.spin_polarized
    calculator.cell.kpoint_mp_grid = parameters.kpoint_mp_grid

    # Prevent the bulk cell from rotating
    calculator.cell.cell_constraints = "0 0 0\n0 0 0"
    calculator.task = task

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


def _get_constraints(
    atoms: Atoms,
    n_free_layer: int,
) -> list[FixAtoms | FixedLine]:
    heights = _get_heights_sorted(atoms)
    sorted_indices = np.argsort(heights)

    free_layer_indices = [i for i in sorted_indices if i < n_free_layer]
    fixed_layer_indices = [i for i in sorted_indices if i >= n_free_layer]
    return [
        FixAtoms(fixed_layer_indices),
        FixedLine(free_layer_indices, [0, 0, 1]),  # type: ignore bad lib
    ]


def get_surface(
    atom: Atoms,
    slab_direction: tuple[int, int, int],
    *,
    n_fixed_layer: int = 0,
    n_free_layer: int = 0,
    n_vacuum_layer: int = 0,
) -> Atoms:
    n_layer = n_free_layer + n_fixed_layer
    atoms = cast(Atoms, surface(atom, slab_direction, n_layer, 0))
    height_per_layer = _get_height_per_layer(atom, slab_direction)
    add_vacuum(atoms, n_vacuum_layer * height_per_layer)

    atoms.set_constraint(  # type: ignore bad lib
        _get_constraints(atoms, n_free_layer),
    )

    return atoms


def _get_n_vacuum_layers_from_slab(atom: Atoms) -> float:
    cell_height = atom.get_cell().array[-1][-1]

    slab_height = get_n_th_slab_layer_height(atom, 1)
    height_per_layer = get_height_per_slab_layer(atom)

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
    n_slab_layer = get_n_slab_layer(atom)
    atom_height = get_n_th_slab_layer_height(atom, n_th_layer)
    height_per_slab_layer = get_height_per_slab_layer(atom)
    atom_height_before_relax = (n_slab_layer - (n_th_layer - 1)) * height_per_slab_layer
    return atom_height_before_relax - atom_height


def plot_displacement_against_n_free_layer(
    calculators: list[Castep],
    layer_idx: int,
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    """Plot displacement of free atoms from it's initial position."""
    displacement = list[float]()
    n_free_layers = list[int]()
    for calculator in calculators:
        atom = get_calculator_atom(calculator)

        displacement.append(
            get_atom_displacement(atom, layer_idx),
        )
        fixed_layers = set[int]()
        for constraint in atom.constraints:  # type: ignore bad lib
            if isinstance(constraint, FixAtoms):
                fixed_layers.update(constraint.index)  # type: ignore bad lib
        n_fixed_layer = len(fixed_layers)
        n_tot_layer = get_n_slab_layer(atom)
        n_free_layers.append(n_tot_layer - n_fixed_layer)

    fig, ax, line = plot_data_comparison(
        ("N Free Layers", np.array(n_free_layers), None),
        ("Displacement", np.array(displacement), "m"),
        ax=ax,
    )
    return fig, ax, line
