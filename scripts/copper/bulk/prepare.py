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

if __name__ == "__main__":
    bulk_copper = bulk("Cu", "fcc", 3.8)  #
    data_path = Path("data")

    calculators = list[Castep]()

    for n_k_points in range(10):
        config = CastepConfig(data_path / f"bulk_{n_k_points}", "bulk")
        params = BulkOptimizationParams(n_k_points=n_k_points)
        calculator = get_bulk_optimization_calculator(bulk_copper, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)

    prepare_submit_all_script(calculators, data_path)
    copy_files_to_hpc(data_path, PosixPath("copper/bulk"))
