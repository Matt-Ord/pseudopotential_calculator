from __future__ import annotations

from pathlib import Path, PosixPath
from typing import TYPE_CHECKING

from ase.build import bulk  # type: ignore bad library

from pseudopotential_calculator.calculations.bulk import (
    BulkOptimizationParams,
    XCFunctional,
    get_bulk_optimization_calculator,
)
from pseudopotential_calculator.castep import Castep, CastepConfig
from pseudopotential_calculator.hpc import (
    prepare_calculator_with_submit_script,
    prepare_submit_all_script,
)
from pseudopotential_calculator.scripting import maybe_copy_files_to_hpc
from pseudopotential_calculator.util import (
    plot_atoms,
    prepare_clean_directory,
    save_fig,
)

if TYPE_CHECKING:
    from ase import Atoms


def _prepare_k_points_convergence(
    atom: Atoms,
    data_path: Path,
    *,
    xc_functional: XCFunctional = "PBE",
    spin_polarized: bool = True,
) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for n_k_points in range(1, 11):
        config = CastepConfig(data_path / f"bulk_{n_k_points}", "bulk")
        params = BulkOptimizationParams(
            n_k_points=n_k_points,
            xc_functional=xc_functional,
            spin_polarized=spin_polarized,
        )
        calculator = get_bulk_optimization_calculator(atom, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)

    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


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
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


SAVE_DIR = Path(__file__).parent / "figures"


def _visualize_initial_atoms(atom: Atoms) -> None:
    fig, _ = plot_atoms(atom, radii=0.3, rotation=(10, 0, 0))
    plot_name = "initial_arrangement.png"
    save_fig(fig, SAVE_DIR / plot_name)


K_POINTS_PATH_PBE = PosixPath("data/copper/bulk/k_points_PBE")
K_POINTS_PATH_WC = PosixPath("data/copper/bulk/k_points_WC")
K_POINTS_PATH_SP = PosixPath("data/copper/bulk/k_points_SP")
ENERGY_CUTOFF_PATH = PosixPath("data/copper/bulk/cutoff_energy")

if __name__ == "__main__":
    bulk_copper = bulk("Cu", "fcc", 3.8)

    _prepare_k_points_convergence(bulk_copper, K_POINTS_PATH_WC, xc_functional="WC")
    _prepare_k_points_convergence(bulk_copper, K_POINTS_PATH_SP, spin_polarized=True)
    _prepare_k_points_convergence(bulk_copper, K_POINTS_PATH_PBE, xc_functional="PBE")
    _prepare_cutoff_energy_convergence(bulk_copper, ENERGY_CUTOFF_PATH)

    _visualize_initial_atoms(bulk_copper)
