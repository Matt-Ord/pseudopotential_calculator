from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Self, cast

import matplotlib.pyplot as plt
import numpy as np
from ase import Atom, Atoms
from matplotlib import cm

from pseudopotential_calculator.calculations.generic import (
    OptimizationParamsBase,
    repeat_cell,
)
from pseudopotential_calculator.calculations.slab import get_n_th_slab_layer_height
from pseudopotential_calculator.castep import (
    CastepConfig,
    get_calculator_atom,
    get_default_calculator,
)

if TYPE_CHECKING:
    from ase.calculators.castep import Castep
    from matplotlib.figure import Figure


class AdsorbateOptimizationParams(OptimizationParamsBase):
    """Parameters of a adsorbate calculation."""

    @property
    def kpoint_mp_grid(self: Self) -> str:
        return f"{self.n_k_points} {self.n_k_points} {1}"


def get_adsorbate_optimization_calculator(
    atoms: Atoms,
    parameters: AdsorbateOptimizationParams,
    config: CastepConfig,
    task: Literal["Energy", "GeometryOptimization"] = "Energy",
) -> Castep:
    calculator = get_default_calculator(config)

    calculator.param.xc_functional = parameters.xc_functional
    calculator.param.cut_off_energy = parameters.cut_off_energy
    calculator.param.spin_polarized = parameters.spin_polarized
    calculator.cell.kpoint_mp_grid = parameters.kpoint_mp_grid

    # Prevent the adsorbate cell from rotating
    calculator.cell.cell_constraints = "1 1 1\n0 0 0"
    calculator.task = task

    calculator.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculator
    return calculator


def get_top_position(atoms: Atoms) -> list[float]:
    positions = atoms.positions  # type: ignore bad lib
    return max(positions, key=lambda x: x[2])  # type: ignore bad lib


def get_top_layer_z_vector(atoms: Atoms) -> list[float]:
    z = get_n_th_slab_layer_height(atoms, 0)
    return [0, 0, z]


def prepare_na_atom(
    slab_position: list[float],
    rel_position: list[float] | None = None,
) -> Atom:
    if rel_position is None:
        rel_position = [0, 0, 1.0]
    na_position = [slab + rel for slab, rel in zip(slab_position, rel_position)]
    return Atom("Na", position=na_position)


def add_atom_onto_slab(atom: Atom, slab: Atoms) -> Atoms:
    slab.append(atom)  # type: ignore bad lib
    return slab


def prepare_adsorbate(atom: Atom, slab: Atoms, width: int) -> Atoms:
    repeated_slab = repeat_cell(slab, (width, width, 1))
    add_atom_onto_slab(atom, repeated_slab)
    return repeated_slab


def get_basis_vectors(slab: Atoms, n: int = 6) -> list[list[float]]:
    lattice_parameters = cast(list[list[float]], slab.cell)
    x = lattice_parameters[0]
    y = lattice_parameters[1]
    return [[i / n for i in x], [i / n for i in y]]


def plot_energy_against_position(
    calculators: list[Castep],
) -> Figure:
    """Plot displacement of free atoms from it's initial position."""
    positions = list[list[float]]()
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
