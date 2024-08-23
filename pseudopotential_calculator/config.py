# config.py

# Lattice parameter for Copper
from dataclasses import dataclass
from pathlib import Path

# CASTEP settings


@dataclass
class CASTEPConfig:
    # k-point grid
    k_range = range(2, 10)
    k_point = 6
    k_points: int = 6, 6, 6
    # atom
    atom_name: str
    lattice_parameter: float = 3.8  # Angstroms
    lattice_type: str = "fcc"
    # calculation
    cut_off_energy: float = 340  # eV
    xc_functional: str = "PBE"
    castep_command: str = "castep.serial"
    # save data
    calculation_directory: Path
    file_name: str
