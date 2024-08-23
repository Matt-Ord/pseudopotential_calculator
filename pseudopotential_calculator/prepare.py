from ast import Dict
from pathlib import Path

from ase.build import bulk
from ase.calculators.castep import Castep
from numpy import nonzero

from pseudopotential_calculator.config import CASTEPConfig

def get_castep_calculator(

    directory: str,
    label: str,
    k: int,
    config: CASTEPConfig,
) -> None:
    """Set up the CASTEP calculation parameters using the values from config."""
    calculation = Castep(
        directory=directory,
        label=label,
        castep_command=config.castep_command,  # Use the command from config
        check_castep_version=True,
        keyword_tolerance=0,
    )

    # Use the configuration values
    calculation.param.xc_functional = (
        config.xc_functional
    )  # Use the XC functional from config
    calculation.cut_off_energy = (
        config.cut_off_energy
    )  # Use the cutoff energy from config
    calculation.cell.kpoint_mp_grid = f"{k} {k} {k}"  # Use k-points
    calculation.cell.cell_constraints = "1 1 1\n0 0 0"  # type: ignore
    calculation.task = "GeometryOptimization"  # type: ignore

    # Optional: set additional settings if needed
    calculation._track_output = False  # type: ignore # noqa: SLF001

# TODO delete submit_script files and use environmental variables
def prepare_submit_script(
    submit_script_from_path: str = "scripts/submit_script.sh" ,
    config: CASTEPConfig,
) -> None:

    """Read the submit script, replace placeholders, and write it to the target location."""
    with open(submit_script_from_path) as file:
            content = file.read()
    for k in config.k_range:
        directory = f"data/bulk_cu_{k}x{k}x{k}"
        submit_script_to_path = os.path.join(directory, "submit_job.sh")
        replacements = {"bulk_cu_nxnxn": f"bulk_cu_{k}x{k}x{k}"}
        for n, k in replacements.items():
            content = content.replace(n, k)

        with open(submit_script_to_path, "w") as f:
            f.write(content)

    """Write the submit_all_script into the data directory."""
    submit_all_script_to_path = os.path.join("data", "submit_all.sh")
    submit_all_script_from_path = os.path.join("scripts/submit_all_script.sh" )

    with open(submit_all_script_from_path, "r") as template_file:
        content = template_file.read()

    with open(submit_all_script_to_path, "w") as f:

        f.write(content)

def prepare_calculations(atom, config:CASTEPConfig) -> None:
    """Prepare and run calculations for the given range of k."""
    for k in config.k_range:
        directory = f"data/bulk_cu_{k}x{k}x{k}"
        label = f"bulk_cu_{k}x{k}x{k}"
        calculation = get_castep_calculator(directory, label, k)
        calculation.set_atoms(atom)
        atom.set_calculator(calculation)

        if calculation.dryrun_ok():
            print("dry run ok")
            calculation.prepare_input_files()
        else:
            msg = "dryrun failed"
            raise Exception(msg)
