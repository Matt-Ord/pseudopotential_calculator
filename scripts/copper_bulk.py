from ase.build import bulk
from ase.calculators.castep import Castep

# Define bulk copper structure
bulk_copper = bulk("Cu", "fcc", a=3.615)  # FCC copper, a = 3.615 Å

# Set up the CASTEP calculator for bulk copper
directory = "data/bulk_cu"
calculation = Castep(
    directory=directory,
    label="bulk_cu",
    castep_command="castep.serial",
    check_castep_version=True,
)
calculation.xc_functional = "PBE"
calculation.cut_off_energy = 340  # in eV
calculation.cell.kpoint_mp_grid = "1 1 1"
calculation.task = "GeometryOptimization"
calculation._track_output = True
calculation._try_reuse = True
calculation._pedantic = True
calculation._rename_existing_dir = False
calculation.param.reuse = True
calculation.param.num_dump_cycles = 0
calculation._export_settings = True

# Attach the calculator to the atoms object
calculation.set_atoms(bulk_copper)


if calculation.dryrun_ok():
    print(f"{calculation._label} : {calculation.get_potential_energy()} ")
else:
    print("Found error in input")
    print(calculation._error)
