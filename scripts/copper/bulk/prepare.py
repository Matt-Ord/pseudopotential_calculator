from pathlib import Path, PosixPath

from ase.build import bulk  # type: ignore bad library

from pseudopotential_calculator.calculations.bulk import (
    BulkOptimizationParams,
    get_bulk_optimization_calculator,
)
from pseudopotential_calculator.castep import Castep, CastepConfig
from pseudopotential_calculator.hpc import (
    copy_files_to_hpc,
    prepare_calculator_with_submit_script,
    prepare_submit_all_script,
)
from pseudopotential_calculator.util import prepare_clean_directory

if __name__ == "__main__":
    bulk_copper = bulk("Cu", "fcc", 3.8)  #
    data_path_k = Path("data/copper/bulk_k")
    data_path_ce = Path("data/copper/bulk_ce")
    prepare_clean_directory(data_path_k)

    prepare_clean_directory(data_path_ce)

    calculators = list[Castep]()

    # k points convergence test
    for n_k_points in range(1, 11):
        config = CastepConfig(data_path_k / f"bulk_{n_k_points}", "bulk")
        params = BulkOptimizationParams(n_k_points=n_k_points)
        calculator = get_bulk_optimization_calculator(bulk_copper, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    # cut off energy convergence test
    for ce in range(300, 750, 50):
        config = CastepConfig(data_path_ce / f"bulk_{ce}", "bulk")
        params = BulkOptimizationParams(n_k_points=11, cut_off_energy=ce)
        calculator = get_bulk_optimization_calculator(bulk_copper, params, config)
        prepare_calculator_with_submit_script(calculator)

    prepare_submit_all_script(calculators, data_path_k)

    prepare_submit_all_script(calculators, data_path_ce)
    copy_files_to_hpc(data_path_k, PosixPath("copper/bulk"))

    copy_files_to_hpc(data_path_ce, PosixPath("copper/bulk"))
