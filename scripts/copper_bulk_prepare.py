import os

from ase.build import bulk
from ase.calculators.castep import Castep

# Define bulk copper structure
bulk_copper = bulk("Cu", "fcc", a=3.8)  # FCC copper, a = 3.615
# Create castep file with submit script to be submitted to hpc for calculation
for k in range(1, 11):
    directory = f"data/bulk_cu_{k}x{k}x{k}"
    label = f"bulk_cu_{k}x{k}x{k}"
    submit_script_to_path = os.path.join(directory, "submit_job.sh")
    submit_script_from_path = "scripts/submit_script.sh"
    replacements = {
        "bulk_cu_nxnxn": f"bulk_cu_{k}x{k}x{k}",
    }
    calculation = Castep(
        directory=directory,
        label=label,
        castep_command="castep.serial",
        check_castep_version=True,
        keyword_tolerance=0,
    )
    calculation.xc_functional = "PBE"
    calculation.cut_off_energy = 340  # in eV
    calculation.cell.kpoint_mp_grid = f"{k} {k} {k}"
    # Controlling the angles within the cell so that it's FCC all the time during simulation
    calculation.cell.cell_constraints = "1 1 1\n0 0 0"
    calculation.task = "GeometryOptimization"

    calculation._track_output = False  # type: ignore : bad library type
    calculation._try_reuse = True
    calculation._pedantic = True
    calculation._rename_existing_dir = False
    calculation.param.reuse = True
    calculation.param.num_dump_cycles = 0
    calculation._export_settings = True
    # Attach the calculator to the atoms object
    calculation.set_atoms(bulk_copper)
    bulk_copper.set_calculator(calculation)

    if calculation.dryrun_ok():
        print("ok")

        calculation.prepare_input_files()
        # Write the submit script

        with open(submit_script_from_path) as file:
            content = file.read()

        for n, k in replacements.items():
            content = content.replace(n, k)

        with open(submit_script_to_path, "w") as f:
            f.write(content)
    else:
        msg = "dryrun failed"
        raise Exception(msg)

    # write a submit_all_script into data file
    submit_all_script_from_path = "scripts/submit_all_script.sh"
    submit_all_script_to_path = os.path.join("data", "submit_all.sh")
    with open(submit_all_script_from_path) as file:
        content = file.read()

    with open(submit_all_script_to_path, "w") as f:
        f.write(content)


import os

from ase.build import bulk
from ase.calculators.castep import Castep


def create_bulk_copper_structure():
    """Create and return the bulk copper structure."""
    return bulk("Cu", "fcc", a=3.8)


def setup_calculation(directory, label, k):
    """Set up the CASTEP calculation parameters."""
    calculation = Castep(
        directory=directory,
        label=label,
        castep_command="castep.serial",
        check_castep_version=True,
        keyword_tolerance=0,
    )
    calculation.xc_functional = "PBE"
    calculation.cut_off_energy = 340  # in eV
    calculation.cell.kpoint_mp_grid = f"{k} {k} {k}"
    calculation.cell.cell_constraints = "1 1 1\n0 0 0"
    calculation.task = "GeometryOptimization"

    calculation._track_output = False  # type: ignore : bad library type
    calculation._try_reuse = True
    calculation._pedantic = True
    calculation._rename_existing_dir = False
    calculation.param.reuse = True
    calculation.param.num_dump_cycles = 0
    calculation._export_settings = True

    return calculation


def write_submit_script(
    submit_script_from_path,
    submit_script_to_path,
    replacements,
) -> None:
    """Read the submit script, replace placeholders, and write it to the target location."""
    with open(submit_script_from_path) as file:
        content = file.read()

    for n, k in replacements.items():
        content = content.replace(n, k)

    with open(submit_script_to_path, "w") as f:
        f.write(content)


def prepare_calculations(k_range) -> None:
    """Prepare and run calculations for the given range of k."""
    bulk_copper = create_bulk_copper_structure()

    for k in k_range:
        directory = f"data/bulk_cu_{k}x{k}x{k}"
        label = f"bulk_cu_{k}x{k}x{k}"
        submit_script_to_path = os.path.join(directory, "submit_job.sh")
        submit_script_from_path = "scripts/submit_script.sh"
        replacements = {"bulk_cu_nxnxn": f"bulk_cu_{k}x{k}x{k}"}

        calculation = setup_calculation(directory, label, k)
        calculation.set_atoms(bulk_copper)
        bulk_copper.set_calculator(calculation)

        if calculation.dryrun_ok():
            print("ok")
            calculation.prepare_input_files()
            write_submit_script(
                submit_script_from_path,
                submit_script_to_path,
                replacements,
            )
        else:
            msg = "dryrun failed"
            raise Exception(msg)


def write_submit_all_script() -> None:
    """Write the submit_all_script into the data directory."""
    submit_all_script_from_path = "scripts/submit_all_script.sh"
    submit_all_script_to_path = os.path.join("data", "submit_all.sh")

    with open(submit_all_script_from_path) as file:
        content = file.read()

    with open(submit_all_script_to_path, "w") as f:
        f.write(content)


if __name__ == "__main__":
    k_range = range(1, 11)
    prepare_calculations(k_range)
    write_submit_all_script()
