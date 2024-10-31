from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

from pseudopotential_calculator.atoms import (
    Vector,
    add_vectors,
    append_atom_to_atoms,
    get_basis_vectors,
    repeat_cell,
)
from pseudopotential_calculator.calculations.slab import get_height_of_layer
from pseudopotential_calculator.castep import (
    get_calculator_atom,
)

if TYPE_CHECKING:
    from ase import Atom, Atoms
    from ase.calculators.castep import Castep
    from matplotlib.figure import Figure


def get_top_position(atoms: Atoms) -> tuple[float, float, float]:
    positions = atoms.positions  # type: ignore bad lib
    return tuple(max(positions, key=lambda x: x[2]))  # type: ignore bad lib


def _get_top_layer_z_vector(atoms: Atoms) -> tuple[float, float, float]:
    z = get_height_of_layer(atoms, 0)
    return (0, 0, z)


# TODO add vacuum
# do relatie till the end
def prepare_slab_with_adsorbate(atom: Atom, slab: Atoms, width: int) -> Atoms:
    repeated_slab = repeat_cell(slab, (width, width, 1))
    return append_atom_to_atoms(atom, repeated_slab)


def _get_scaled_basis_vectors(
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


def prepare_adsorbate_positions(
    slab: Atoms,
    heights_above_surface: list[float],
    *,
    scale_factor: int,
) -> list[tuple[float, float, float]]:
    adsorbate_positions = list[tuple[float, float, float]]()

    basis_vectors = _get_scaled_basis_vectors(slab, scale_factor=scale_factor)
    z = _get_top_layer_z_vector(slab)
    x = basis_vectors[0]
    y = basis_vectors[1]

    i_vals, j_vals, k_vals = np.meshgrid(
        np.arange(scale_factor + 1),
        np.arange(scale_factor + 1),
        heights_above_surface,
        indexing="ij",
    )

    i_flat = i_vals.ravel()
    j_flat = j_vals.ravel()
    k_flat = k_vals.ravel()

    for idx in range(len(i_flat)):
        i, j, k = i_flat[idx], j_flat[idx], k_flat[idx]

        if i >= j:
            slab_position = (
                i * x[0] + j * y[0] + z[0],
                i * x[1] + j * y[1] + z[1],
                i * x[2] + j * y[2] + z[2],
            )
            adsorbate_position = add_vectors(slab_position, (0, 0, k))
            adsorbate_positions.append(adsorbate_position)
    return adsorbate_positions
