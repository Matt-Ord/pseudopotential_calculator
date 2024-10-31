from pathlib import Path, PosixPath

import numpy as np
from ase import Atoms

from pseudopotential_calculator.atoms import repeat_cell
from pseudopotential_calculator.calculations.slab import (
    get_surface,
)
from pseudopotential_calculator.castep import (
    Castep,
    CastepConfig,
    CastepParams,
    get_calculator,
    load_calculator_atoms,
)
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

VACUUM_HEIGHT_PATH = Path("data/copper/slab/vacuum_height")
FREE_LAYER_PATH = Path("data/copper/slab/free_layer")
SLAB_WITH_PATH = Path("data/copper/slab/width")


def _prepare_vacuum_height_convergence(atoms: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for vacuum_height in range(1, 7):
        slab = get_surface(
            atoms,
            (1, 1, 1),
            n_fixed_layer=5,
            vacuum_height=vacuum_height,
        )

        config = CastepConfig(
            data_path / f"{vacuum_height}",
            f"s_{vacuum_height}_vacuum_height",
        )
        params = CastepParams(n_k_points=(12, 12, 1), xc_functional="WC")
        calculator = get_calculator(slab, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


def _prepare_free_layer_convergence(atoms: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)
    for n_free_layer in range(1, 6):
        n_fixed_layer = n_free_layer + 1
        slab = get_surface(
            atoms,
            (1, 1, 1),
            n_fixed_layer=n_fixed_layer,
            n_free_layer=n_free_layer,
            vacuum_height=6,
        )

        config = CastepConfig(
            data_path / f"{n_free_layer}",
            f"s_{n_free_layer}_free_layer",
        )
        params = CastepParams(
            n_k_points=(12, 12, 1),
            xc_functional="WC",
        )
        calculator = get_calculator(
            slab,
            params,
            config,
        )
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


def _prepare_n_widths_slab(
    atoms: Atoms,
    data_path: Path,
) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    slab = get_surface(
        atoms,
        (1, 1, 1),
        n_fixed_layer=4,
        n_free_layer=3,
        vacuum_height=10,
    )

    for width in range(1, 8):
        n_width_slab = repeat_cell(slab, (width, width, 1))
        config = CastepConfig(
            data_path / f"{width}",
            f"s_{width}_width",
        )
        params = CastepParams(
            n_k_points=(int(np.ceil(12 / width)), int(np.ceil(12 / width)), 1),
            xc_functional="WC",
            max_scf_cycles=100,
            cut_off_energy=600,
        )
        calculator = get_calculator(n_width_slab, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


if __name__ == "__main__":
    bulk_config = CastepConfig(Path("data/copper/bulk/k_points_WC/bulk_10"), "bulk")
    bulk_copper = load_calculator_atoms(bulk_config)

    _prepare_vacuum_height_convergence(bulk_copper, VACUUM_HEIGHT_PATH)

    _prepare_free_layer_convergence(bulk_copper, FREE_LAYER_PATH)

    _prepare_n_widths_slab(bulk_copper, SLAB_WITH_PATH)

    slab = get_surface(
        bulk_copper,
        (1, 1, 1),
        n_fixed_layer=5,
        vacuum_height=8,
    )

    fig, _ = plot_atoms(slab, radii=0.3, rotation=(-90, +30, 0))
    SAVE_DIR = Path(__file__).parent / "figures"
    plot_name = "initial_arrangement.png"
    save_fig(fig, SAVE_DIR / plot_name)
