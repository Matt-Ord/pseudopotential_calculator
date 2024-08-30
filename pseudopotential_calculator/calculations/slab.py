from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Self, cast

from ase import Atoms
from ase.build import (
    surface,  # type: ignore  # noqa: PGH003
)
from ase.constraints import FixAtoms

from pseudopotential_calculator.calculations.generic import OptimizationParamsBase
from pseudopotential_calculator.castep import CastepConfig, get_default_calculator

if TYPE_CHECKING:
    from ase.calculators.castep import Castep


@dataclass
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
    thickness: float,
    vacuum: float,
) -> Atoms:
    # Call the surface function and cast the result to Atoms
    slab = surface(atom, slab_direction, thickness, vacuum=vacuum)  # type: ignore bad library
    return cast(Atoms, slab)
