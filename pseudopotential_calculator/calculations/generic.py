from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Self

if TYPE_CHECKING:
    from ase import Atoms


def repeat_cell(cell: Atoms, amount: tuple[int, int, int]) -> Atoms:
    return cell.repeat(amount)  # type: ignore bad lib


@dataclass
class OptimizationParamsBase(ABC):
    """Base class for optimization parameters."""

    n_k_points: int = field(default=1, kw_only=True)
    cut_off_energy: float = field(default=600, kw_only=True)
    xc_functional: str = field(
        default="PBE",
        kw_only=True,
    )  # Adjusted to str for generalization
    spin_polarized: bool = field(default=False, kw_only=True)
    max_scf_cycles: int = 30

    @property
    @abstractmethod
    def kpoint_mp_grid(self: Self) -> str:
        """Method to return a string of k points to be implemented by subclasses."""
