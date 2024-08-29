from pathlib import Path, PosixPath

from ase import Atoms
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


def _prepare_k_points_convergence(atom: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for n_k_points in range(1, 11):
        config = CastepConfig(data_path / f"bulk_{n_k_points}", "bulk")
        params = BulkOptimizationParams(n_k_points=n_k_points)
        calculator = get_bulk_optimization_calculator(atom, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)

    prepare_submit_all_script(calculators, data_path)
    copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


def _prepare_k_points_convergence_WC(atom: Atoms, data_path: Path) -> None:  # noqa: N802
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for n_k_points in range(1, 11):
        config = CastepConfig(data_path / f"bulk_{n_k_points}", "bulk")
        params = BulkOptimizationParams(
            n_k_points=n_k_points,
            cut_off_energy=600,
            xc_functional="WC",
        )
        calculator = get_bulk_optimization_calculator(atom, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)

    prepare_submit_all_script(calculators, data_path)
    copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


def _prepare_cutoff_energy_convergence(atom: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for cutoff_energy in range(300, 750, 50):
        config = CastepConfig(data_path / f"bulk_{cutoff_energy}", "bulk")
        params = BulkOptimizationParams(n_k_points=11, cut_off_energy=cutoff_energy)
        calculator = get_bulk_optimization_calculator(atom, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)

    prepare_submit_all_script(calculators, data_path)
    copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


K_POINTS_PATH = Path("data_min/copper/bulk/k_points")
K_POINTS_PATH_WC = Path("data_min/copper/bulk/k_points_WC")
ENERGY_CUTOFF_PATH = Path("data_min/copper/bulk/cutoff_energy")

if __name__ == "__main__":
    bulk_copper = bulk("Cu", "fcc", 3.8)

    _prepare_k_points_convergence_WC(bulk_copper, K_POINTS_PATH_WC)
    _prepare_cutoff_energy_convergence(bulk_copper, ENERGY_CUTOFF_PATH)
    _prepare_k_points_convergence(bulk_copper, K_POINTS_PATH)
