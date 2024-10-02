from __future__ import annotations

from pathlib import PosixPath
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from ase import Atoms
from ase.build import (
    add_vacuum,  # type: ignore bad library
    surface,  # type: ignore bad library
)
from ase.constraints import FixAtoms, FixedLine

from pseudopotential_calculator.castep import (
    get_atom_potential_energy,
    get_calculator_atom,
    get_calculator_old_atoms,
    load_all_calculators,
)
from pseudopotential_calculator.util import plot_data_comparison

if TYPE_CHECKING:
    from ase.calculators.castep import Castep
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

PLAIN_SLAB_PATH = PosixPath("data/copper/slab/width")


def _get_unsorted_heights(atoms: Atoms) -> np.ndarray[Any, np.dtype[np.float64]]:
    heights = cast(np.ndarray[Any, np.dtype[np.float64]], atoms.positions[:, 2])  # type: ignore bad lib
    return np.unique(heights)


def _get_unsorted_widths(atoms: Atoms) -> np.ndarray[Any, np.dtype[np.float64]]:
    widths = cast(np.ndarray[Any, np.dtype[np.float64]], atoms.positions[:, 1])  # type: ignore bad lib
    return np.unique(widths)


def _get_sorted_heights(atoms: Atoms) -> np.ndarray[Any, np.dtype[np.float64]]:
    return _get_unsorted_heights(atoms)[::-1]


def get_height_of_layer(atoms: Atoms, idx: int) -> float:
    heights = _get_sorted_heights(atoms)
    return heights[idx]


def get_height_per_layer(atom: Atoms) -> float:
    heights = _get_sorted_heights(atom)
    return heights[-2] - heights[-1]


def get_n_slab_layer(atoms: Atoms) -> int:
    heights = _get_unsorted_heights(atoms)
    return len(heights)


def get_slab_width(atoms: Atoms) -> int:
    # n_slab_atoms / n_atoms_per_repeat gives the number of slab units set to be repeated.
    # take square root to get width since the repeat s set to be the same in x and y directions
    width = int(len(_get_unsorted_widths(atoms)) / 3)
    return width  # noqa: RET504


def _get_vacuum_height_from_slab(atoms: Atoms) -> float:
    cell_height = atoms.get_cell().array[-1][-1]

    slab_height = get_height_of_layer(atoms, 0)

    return cell_height - slab_height


def _get_sorted_calculators(
    calculators: list[Castep],
    n_adsorbate_atom: int,
) -> list[Castep]:
    return sorted(
        calculators,
        key=lambda slab_calculator: get_slab_width(
            get_calculator_atom(slab_calculator),
        ),
    )


def plot_energy_against_slab_width(
    slab_calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    widths = list[int]()
    energies = list[float]()

    plain_slab_calculators = load_all_calculators(PLAIN_SLAB_PATH)
    sorted_plain_slab_calculators = _get_sorted_calculators(
        plain_slab_calculators,
        n_adsorbate_atom=0,
    )

    sorted_slab_calculators = _get_sorted_calculators(
        slab_calculators,
        n_adsorbate_atom=1,
    )

    for slab_calculator, plain_slab_calculator in zip(
        sorted_slab_calculators,
        sorted_plain_slab_calculators,
    ):
        slab = get_calculator_atom(slab_calculator)
        plain_slab = get_calculator_atom(plain_slab_calculator)
        width = get_slab_width(slab)
        energy = get_atom_potential_energy(slab) - get_atom_potential_energy(plain_slab)
        energies.append(energy)
        widths.append(width)

    widths = np.array(widths)
    energies = np.array(energies)

    fig, ax, original_line = plot_data_comparison(
        ("Width", widths, "atoms"),
        ("Energy", energies, "eV"),
        ax=ax,
    )

    return fig, ax, original_line


def _get_constraints(
    atoms: Atoms,
    n_free_layer: int,
) -> list[FixAtoms | FixedLine]:
    heights = _get_unsorted_heights(atoms)
    indices = np.argsort(heights)[::-1]

    free_layer_indices = indices[:n_free_layer]
    fixed_layer_indices = indices[n_free_layer:]
    return [
        FixAtoms(fixed_layer_indices),
        FixedLine(free_layer_indices, [0, 0, 1]),  # type: ignore bad lib
    ]


def get_surface(
    atoms: Atoms,
    slab_direction: tuple[int, int, int],
    *,
    n_fixed_layer: int = 0,
    n_free_layer: int = 0,
    vacuum_height: float = 1,  # unit: Angstrom
) -> Atoms:
    n_layer = n_free_layer + n_fixed_layer
    slab = cast(Atoms, surface(atoms, slab_direction, n_layer, 0))
    add_vacuum(
        slab,
        vacuum_height,
    )

    slab.set_constraint(  # type: ignore bad lib
        _get_constraints(slab, n_free_layer),
    )

    return slab


def plot_energy_against_vacuum_height(
    calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    """Plot energy against vacuum height.

    This assumes the vacuum is layered up in z direction.

    Also assumes the bottom of the slab is at the bottom of the cell.
    """
    energies = list[float]()
    vacuum_heights = list[float]()
    for calculator in calculators:
        atom = get_calculator_atom(calculator)
        old_atom = get_calculator_old_atoms(calculator)
        energies.append(get_atom_potential_energy(atom))
        vacuum_height = _get_vacuum_height_from_slab(old_atom)
        vacuum_heights.append(vacuum_height)

    return plot_data_comparison(
        ("Vacuum Height", np.array(vacuum_heights), "Angstrom"),
        ("Energy", np.array(energies), "eV"),
        ax=ax,
    )


def get_displacement_of_layer(
    idx: int,
    *,
    atoms: Atoms,
    initial_atoms: Atoms,
) -> float:
    n_slab_layers = get_n_slab_layer(atoms)
    height_per_layer = get_height_per_layer(initial_atoms)
    if idx >= n_slab_layers or idx < -n_slab_layers:
        raise IndexError
    idx = idx % n_slab_layers
    final_height = get_height_of_layer(atoms, idx)
    initial_height = (n_slab_layers - idx - 1) * height_per_layer
    return final_height - initial_height


def plot_displacement_against_n_free_layer(
    calculators: list[Castep],
    *,
    idx_0: int,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    """Plot displacement of free atoms from it's initial position.

    The layer index is defined such that index 0 has the largest z value
    and index -1 has the smallest z value.
    """
    displacements = list[float]()
    n_free_layers = list[int]()

    for calculator in calculators:
        atoms = get_calculator_atom(calculator)
        initial_atoms = get_calculator_old_atoms(calculator)
        n_layer = get_n_slab_layer(atoms)

        idx = idx_0
        if idx >= n_layer or idx < -n_layer:
            continue

        displacement = get_displacement_of_layer(
            idx,
            atoms=atoms,
            initial_atoms=initial_atoms,
        )
        displacements.append(displacement)
        fixed_layers = set[int]()
        for constraint in atoms.constraints:  # type: ignore bad lib
            if isinstance(constraint, FixAtoms):
                fixed_layers.update(constraint.index)  # type: ignore bad lib
        n_fixed_layer = len(fixed_layers)
        n_free_layers.append(n_layer - n_fixed_layer)

    fig, ax, line = plot_data_comparison(
        ("N Free Layers", np.array(n_free_layers), None),
        ("Displacement", np.array(displacements), "Å"),
        ax=ax,
    )
    return fig, ax, line
