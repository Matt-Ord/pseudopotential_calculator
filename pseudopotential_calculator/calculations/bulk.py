from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Self, cast

from ase import Atoms
from surface_potential_analysis.util.plot import get_figure

from pseudopotential_calculator.castep import CastepConfig, get_default_calculator

if TYPE_CHECKING:
    from pathlib import Path

    from ase.calculators.castep import Castep
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D


@dataclass
class BulkOptimizationParams:
    """Parameters of a bulk calculation."""

    n_k_points: int = field(default=1, kw_only=True)
    cut_off_energy: float = field(default=340, kw_only=True)
    xc_functional: str = field(default="PBE", kw_only=True)

    @property
    def kpoint_mp_grid(self: Self) -> str:
        return f"{self.n_k_points} {self.n_k_points} {self.n_k_points}"


def get_bulk_optimization_calculator(
    atoms: Atoms,
    parameters: BulkOptimizationParams,
    config: CastepConfig,
) -> Castep:
    calculation = get_default_calculator(config)

    calculation.param.xc_functional = parameters.xc_functional
    calculation.param.cut_off_energy = parameters.cut_off_energy

    calculation.cell.kpoint_mp_grid = parameters.kpoint_mp_grid
    # Prevent the bulk cell from rotating
    calculation.cell.cell_constraints = "1 1 1\n0 0 0"
    calculation.task = "GeometryOptimization"

    calculation.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculation
    return calculation


def bond_length_convergence_test(
    directory: Path,
) -> tuple[list[int], list[list[float]]]:
    k_values = []
    bond_lengths = []
    previous_bond_length = None
    is_converged = False

    for n in CASTEPConfig.k_range:
        atoms = read(f"./data/bulk_cu_{n}x{n}x{n}/bulk_cu_{n}x{n}x{n}.castep")
        # Get cell parameters using deprecated method
        lengths_and_angles = atoms.get_cell_lengths_and_angles()

        k_values.append(n)
        bond_length = lengths_and_angles[0]
        print(bond_length)
        bond_lengths.append(bond_length)

        if previous_bond_length is not None:
            # Compute the difference in cell parameters
            length_diff = abs(bond_length - previous_bond_length)

            print(f"k = {n}, Bond length = {bond_length}")
            print(f"ΔLengths = {length_diff}")

            if length_diff < convergence_threshold:
                is_converged = True
                print(f"Converged at k = {n}")
                break
        else:
            print(f"k = {n}, Bond length = {bond_length}")
        previous_bond_length = bond_length

    if not is_converged:
        print("The cell parameters did not converge within the tested k-point grids.")

    return k_values, bond_lengths


def _get_unit_cell_displacements_from_calculator(calculator: Castep) -> tuple[float]:
    return cast(Atoms, calculator.atoms).get_cell_lengths_and_angles()[0:3]  # type: ignore unknown


# Function to plot bond length vs k-points
def _plot_bond_length_against_n_k_points(
    n_k_points: list[float],
    bond_lengths: list[float],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    fig, ax = get_figure(ax)

    (line,) = ax.plot(n_k_points, bond_lengths)
    ax.set_xlabel("k-point grid size (n)")
    ax.set_ylabel("Bond Length (Angstrom)")
    ax.set_title("Convergence Test: Bond Length vs k-point grid")
    return fig, ax, line


def plot_bond_length_convergence(
    calculators: list[Castep],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    n_k_points = list[float]()
    bond_lengths = list[float]()
    for calculator in calculators:
        bond_lengths.append(
            _get_unit_cell_displacements_from_calculator(calculator)[0],
        )
        n_k_points.append(calculator.get_kpoints(calculator.atoms)[0])
    return _plot_bond_length_against_n_k_points(n_k_points, bond_lengths, ax=ax)
