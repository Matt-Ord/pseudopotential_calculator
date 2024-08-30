from dataclasses import dataclass, field
from typing import Self


@dataclass
class OptimizationParamsBase:
    """Base class for optimization parameters."""

    n_k_points: int = field(default=1, kw_only=True)
    cut_off_energy: float = field(default=600, kw_only=True)
    xc_functional: str = field(
        default="PBE",
        kw_only=True,
    )  # Adjusted to str for generalization
    spin_polarized: bool = field(default=False, kw_only=True)

    @property
    def kpoint_mp_grid(self: Self) -> str:
        error_message = "This method should be implemented by subclasses."
        raise NotImplementedError(error_message)
