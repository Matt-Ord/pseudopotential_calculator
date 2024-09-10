from __future__ import annotations

from typing import TYPE_CHECKING, Self

from ase import Atom, Atoms

from pseudopotential_calculator.calculations.generic import (
    OptimizationParamsBase,
    repeat_cell,
)
from pseudopotential_calculator.calculations.slab import get_top_position
from pseudopotential_calculator.castep import CastepConfig, get_default_calculator

if TYPE_CHECKING:
    from ase.calculators.castep import Castep


class AdsorbateOptimizationParams(OptimizationParamsBase):
    """Parameters of a adsorbate calculation."""

    @property
    def kpoint_mp_grid(self: Self) -> str:
        return f"{self.n_k_points} {self.n_k_points} {1}"


def get_adsorbate_optimization_calculator(
    atoms: Atoms,
    parameters: AdsorbateOptimizationParams,
    config: CastepConfig,
) -> Castep:
    calculator = get_default_calculator(config)

    calculator.param.xc_functional = parameters.xc_functional
    calculator.param.cut_off_energy = parameters.cut_off_energy
    calculator.param.spin_polarized = parameters.spin_polarized
    calculator.cell.kpoint_mp_grid = parameters.kpoint_mp_grid

    # Prevent the adsorbate cell from rotating
    calculator.cell.cell_constraints = "1 1 1\n0 0 0"
    calculator.task = "GeometryOptimization"

    calculator.set_atoms(atoms)  # type: ignore unknown
    # Temporary fix for bug in ase
    atoms.calc = calculator
    return calculator


def prepare_na_atom(slab: Atoms, rel_position: list[int] | None = None) -> Atom:
    if rel_position is None:
        rel_position = [0, 0, 1]
    rel_position = [0, 0, 1]
    slab_position = get_top_position(slab)
    na_position = [slab + rel for slab, rel in zip(slab_position, rel_position)]
    return Atom("Na", position=na_position)


def add_atom_onto_slab(atom: Atom, slab: Atoms) -> Atoms:
    slab.append(atom)  # type: ignore bad lib
    return slab


def prepare_adsorbate(atom: Atom, slab: Atoms, width: int) -> Atoms:
    repeated_slab = repeat_cell(slab, (width, width, 1))
    add_atom_onto_slab(atom, repeated_slab)
    return repeated_slab
