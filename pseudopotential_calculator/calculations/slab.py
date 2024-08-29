from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Self

from ase.constraints import FixAtoms

from pseudopotential_calculator.castep import CastepConfig, get_default_calculator

if TYPE_CHECKING:
    from ase import Atoms
    from ase.calculators.castep import Castep

    from pseudopotential_calculator.calculations.bulk import XCFunctional


@dataclass
class SlabOptimizationParams:
    """Parameters of a slab calculation."""

    n_k_points: int = field(default=1, kw_only=True)
    cut_off_energy: float = field(default=600, kw_only=True)
    xc_functional: XCFunctional = field(default="PBE", kw_only=True)
    spin_polarized: bool = field(default=False, kw_only=True)

    @property
    def kpoint_mp_grid(self: Self) -> str:
        return f"{self.n_k_points} {self.n_k_points} {1}"


def get_slab_optimization_calculator(
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
