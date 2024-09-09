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


def _get_constraints(
    atoms: Atoms,
    n_free_layer: int,
) -> list[FixAtoms | FixedLine]:
    heights = _get_layer_heights(atoms)
    indices = np.argsort(heights)[::-1]

    free_layer_indices = [i for i in indices if i < n_free_layer]
    fixed_layer_indices = [i for i in indices if i >= n_free_layer]
    return [
        FixAtoms(fixed_layer_indices),
        FixedLine(free_layer_indices, [0, 0, 1]),  # type: ignore bad lib
    ]


def _get_height_per_layer(
    atoms: Atoms,
    slab_direction: tuple[int, int, int],
) -> int:
    # create a slab with two layers,
    # the hight per layer is the difference between the two heights
    slab_single_layer = cast(
        Atoms,
        surface(
            atoms,
            slab_direction,
            layers=2,
            vacuum=0,
        ),
    )

    position_0 = slab_single_layer.positions[0, 2]  # type: ignore unknown
    position_1 = slab_single_layer.positions[1, 2]  # type: ignore unknown
    return np.abs(cast(float, position_0) - cast(float, position_1))


def get_surface(
    atoms: Atoms,
    slab_direction: tuple[int, int, int],
    *,
    n_fixed_layer: int = 0,
    n_free_layer: int = 0,
    n_vacuum_layer: int = 0,
) -> Atoms:
    n_layer = n_free_layer + n_fixed_layer
    slab = cast(Atoms, surface(atoms, slab_direction, n_layer, 0))

    height_per_layer = _get_height_per_layer(atoms, slab_direction)
    add_vacuum(slab, n_vacuum_layer * height_per_layer)

    slab.set_constraint(  # type: ignore bad lib
        _get_constraints(slab, n_free_layer),
    )
    return slab


def _get_heights(atoms: Atoms) -> np.ndarray[Any, np.dtype[np.float64]]:
    return atoms.positions[:, 2]  # type: ignore bad lib


def _get_layer_heights(atoms: Atoms) -> np.ndarray[Any, np.dtype[np.float64]]:
    heights = _get_heights(atoms)
    return np.unique(heights)


def _get_sorted_layer_heights(atoms: Atoms) -> np.ndarray[Any, np.dtype[np.float64]]:
    heights = _get_layer_heights(atoms)
    return np.sort(heights)[::-1]


def _get_height_of_layer(atoms: Atoms, layer_idx: int) -> float:
    heights = _get_sorted_layer_heights(atoms)
    return heights[layer_idx]


def _get_n_vacuum_layers_from_slab(atoms: Atoms) -> float:
    cell_height = atoms.get_cell().array[-1][-1]

    slab_height = _get_height_of_layer(atoms, 1)
    height_per_layer = get_average_height_of_layer(atoms)

    vacuum_height = cell_height - slab_height

    return vacuum_height / height_per_layer


def _get_fixed_atom_indices(atoms: Atoms) -> list[int]:
    fixed_indices = list[int]()
    for constraint in atoms.constraints:  # type: ignore bad lib
        if isinstance(constraint, FixAtoms):
            fixed_indices.append(constraint.index[0])  # type: ignore bad lib
    return fixed_indices


def _get_layer_indices(atoms: Atoms) -> np.ndarray[Any, np.dtype[np.float64]]:
    heights = _get_layer_heights(atoms)
    return cast(np.ndarray[Any, np.dtype[np.float64]], np.argsort(heights))


def get_height_of_fixed_layer(atoms: Atoms) -> float:
    heights = _get_layer_heights(atoms)
    fixed_index = _get_fixed_atom_indices(atoms)
    atom_0 = fixed_index[0]
    atom_1 = fixed_index[1]
    layer_indices = _get_layer_indices(atoms)
    return (heights[atom_0] - heights[atom_1]) / (
        layer_indices[atom_0] - layer_indices[atom_1]
    )


def get_average_height_of_layer(atoms: Atoms) -> float:
    heights = _get_layer_heights(atoms)
    return (max(heights) - min(heights)) / len(heights)


def get_n_slab_layer(atoms: Atoms) -> int:
    return len(_get_layer_heights(atoms))


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


def get_displacement_of_layer(
    atoms: Atoms,
    layer_idx: int,
) -> float:
    n_slab_layers = get_n_slab_layer(atoms)
    height_per_layer = get_height_of_fixed_layer(atoms)
    if layer_idx >= n_slab_layers or layer_idx < -n_slab_layers:
        raise IndexError
    layer_idx = layer_idx % n_slab_layers
    final_height = _get_height_of_layer(atoms, layer_idx)
    initial_height = (n_slab_layers - layer_idx - 1) * height_per_layer

    return final_height - initial_height


def plot_displacement_against_n_free_layer(
    calculators: list[Castep],
    layer_idx_0: int,
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    """Plot displacement of free atoms from it's initial position.

    The layer index is defined such that index 0 has the largest z value
    and index -1 has the smallest z value.
    """
    displacements = list[float]()
    n_free_layers = list[int]()
    n_fixed_layers = list[int]()
    fixed_layers = set[int]()

    for calculator in calculators:
        atom = get_calculator_atom(calculator)
        fixed_layers = set[int]()
        n_layer = get_n_slab_layer(atom)

        layer_idx = layer_idx_0
        if layer_idx >= n_layer or layer_idx < -n_layer:
            continue

        displacement = get_displacement_of_layer(atom, layer_idx)
        displacements.append(displacement)
        for constraint in atom.constraints:  # type: ignore bad lib
            if isinstance(constraint, FixAtoms):
                fixed_layers.update(constraint.index)  # type: ignore bad lib
        n_fixed_layer = len(fixed_layers)
        n_fixed_layers.append(n_fixed_layer)
        n_free_layers.append(n_layer - n_fixed_layer)

    fig, ax, line = plot_data_comparison(
        ("N Free Layers", np.array(n_free_layers), None),
        ("Displacement", np.array(displacements), "Å"),
        ax=ax,
    )
    return fig, ax, line
