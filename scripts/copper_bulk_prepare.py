from ase.build import bulk

from pseudopotential_calculator import prepare
from pseudopotential_calculator.config import CASTEPConfig

config = CASTEPConfig()
bulk_copper = bulk("Cu", "fcc", 3.8)

if __name__ == "__main__":
    prepare.prepare_calculations(bulk_copper, config)
    prepare.prepare_submit_script("scripts/submit_script.sh", config)
