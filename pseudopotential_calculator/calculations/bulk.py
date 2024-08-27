from dataclasses import dataclass, field
from typing import Self

from ase import Atoms
from ase.calculators.castep import Castep

from pseudopotential_calculator.castep import CastepConfig, get_default_calculator


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
