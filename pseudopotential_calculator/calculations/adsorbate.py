from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from ase import Atom, Atoms
from matplotlib import cm

from pseudopotential_calculator.atoms import Vector, get_basis_vectors, repeat_cell
from pseudopotential_calculator.calculations.slab import get_height_of_layer
from pseudopotential_calculator.castep import (
    get_calculator_atom,
)

if TYPE_CHECKING:
    from ase.calculators.castep import Castep
    from matplotlib.figure import Figure


def get_top_position(atoms: Atoms) -> tuple[float, float, float]:
    positions = atoms.positions  # type: ignore bad lib
    return tuple(max(positions, key=lambda x: x[2]))  # type: ignore bad lib


def get_top_layer_z_vector(atoms: Atoms) -> tuple[float, float, float]:
    z = get_height_of_layer(atoms, 0)
    return (0, 0, z)


def get_adsorbate_position(
    slab_position: tuple[float, float, float],
    relative_height: float | None = None,
) -> tuple[float, float, float]:
    if relative_height is None:
        relative_height = 1.0
    return (
        slab_position[0],
        slab_position[1],
        slab_position[2] + relative_height,
    )


def prepare_adsorbate(
    adsorbate: str,
    position: tuple[float, float, float],
) -> Atom:
    return Atom(adsorbate, position)


def get_slab_with_adsorbate(atom: Atom, slab: Atoms) -> Atoms:
    slab.append(atom)  # type: ignore bad lib)
    return slab


def prepare_slab_with_adsorbate(atom: Atom, slab: Atoms, width: int) -> Atoms:
    repeated_slab = repeat_cell(slab, (width, width, 1))
    return get_slab_with_adsorbate(atom, repeated_slab)


def get_scaled_basis_vectors(
    atoms: Atoms,
    *,
    scale_factor: int = 6,
) -> tuple[Vector, Vector, Vector]:
    basis_vectors = get_basis_vectors(atoms)
    x = basis_vectors[0]
    y = basis_vectors[1]
    z = basis_vectors[2]
    return (
        Vector(i / scale_factor for i in x),
        Vector(i / scale_factor for i in y),
        Vector(i / scale_factor for i in z),
    )


def plot_energy_against_position(
    calculators: list[Castep],
) -> Figure:
    """Plot displacement of free atoms from it's initial position."""
    positions = list[tuple[float, float, float]]()
    energies = list[float]()
    for calculator in calculators:
        atom = get_calculator_atom(calculator)
        position = get_top_position(atom)
        positions.append(position)
        energies.append(atom.get_potential_energy())  # type: ignore bad lib
    positions = np.array(positions)
    plt.style.use("_mpl-gallery")

    x = positions[:, 0]
    y = positions[:, 1]
    x, y = np.meshgrid(x, y)
    r = np.sqrt(x**2 + y**2)
    z = np.sin(r)

    # Plot the surface
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})  # type: ignore bad lib
    ax.plot_surface(x, y, z, vmin=z.min() * 2, cmap=cm.Blues)  # type: ignore bad lib

    ax.set(xlabel=[], ylabel=[], zlabel=[])

    return fig
