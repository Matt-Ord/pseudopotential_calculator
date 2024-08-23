import os
from pathlib import Path
from ase.build import bulk
from ase.calculators.castep import Castep

from pseudopotential_calculator.config import CASTEPConfig
from pseudopotential_calculator import prepare

config = CASTEPConfig(
    atom_name= "Cu"
    lattice_parameter= 3.8
    lattice_type= "fcc",
)

atom = bulk("Cu","fcc")



if __name__ == "__main__":



bulk_copper = prepare.create_bulk_structure(config)





    prepare_calculations(config.k_range)  # Use k_range from config
    write_submit_all_script()
