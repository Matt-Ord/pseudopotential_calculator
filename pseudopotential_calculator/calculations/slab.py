from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Self

from pseudopotential_calculator.castep import CastepConfig, get_default_calculator

if TYPE_CHECKING:
    from ase import Atoms
    from ase.calculators.castep import Castep


@dataclass
class SlabOptimizationParams:
    """Parameters of a slab calculation."""

    n_k_points: int = field(default=1, kw_only=True)
    cut_off_energy: float = field(default=600, kw_only=True)
    xc_functional: str = field(default="PBE", kw_only=True)

    @property
    def kpoint_mp_grid(self: Self) -> str:
        return f"{self.n_k_points} {self.n_k_points} {1}"


def get_slab_optimization_calculator(
    atoms: Atoms,
    parameters: SlabOptimizationParams,
    config: CastepConfig,
) -> Castep:
    calculation = get_default_calculator(config)

    calculation.param.xc_functional = parameters.xc_functional
    calculation.param.cut_off_energy = parameters.cut_off_energy
    calculation.param.spinpolarised = "true"
    calculation.param.elec_energy_tol = 1.000000000000000e-06
    calculation.param.geom_energy_tol = 1.000000000000000e-05
    calculation.param.geom_disp_tol = 1.000000000000000e-03
    calculation.cell.kpoint_mp_grid = parameters.kpoint_mp_grid

    # Prevent the bulk cell from rotating
    calculation.cell.cell_constraints = "1 1 1\n0 0 0"
    calculation.task = "GeometryOptimization"

    calculation.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculation
    return calculation
