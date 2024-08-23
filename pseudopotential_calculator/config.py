# config.py

# Lattice parameter for Copper
from dataclasses import dataclass

# CASTEP settings


@dataclass
class CASTEPConfig:
    # k-point grid
    k_range = range(1, 10)
    k_point = 6
    # calculation
    cut_off_energy: float = 340  # eV
    xc_functional: str = "PBE"
    castep_command: str = "castep.serial"
