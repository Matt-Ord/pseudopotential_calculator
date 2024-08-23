from pathlib import Path

from ase.build import bulk  # type: ignore bad library

from pseudopotential_calculator.calculations.bulk import (
    BulkOptimizationParams,
    get_bulk_optimization_calculator,
)
from pseudopotential_calculator.castep import CastepConfig
from pseudopotential_calculator.hpc import prepare_calculator_with_submit_script

if __name__ == "__main__":
    bulk_copper = bulk("Cu", "fcc", 3.8)  #
    data_path = Path("data")

    for n_k_points in range(10):
        config = CastepConfig(data_path / f"bulk_{n_k_points}", "bulk")
        params = BulkOptimizationParams(n_k_points=n_k_points)
        calculator = get_bulk_optimization_calculator(bulk_copper, params, config)
        prepare_calculator_with_submit_script(calculator)
