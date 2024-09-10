from pathlib import Path, PosixPath

from ase import Atoms

from pseudopotential_calculator.calculations.adsorbate import (
    get_top_position,
    prepare_adsorbate,
    prepare_na_atom,
)
from pseudopotential_calculator.calculations.slab import (
    SlabOptimizationParams,
    get_slab_calculator,
    get_surface,
)
from pseudopotential_calculator.castep import (
    Castep,
    CastepConfig,
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

VACUUM_LAYER_PATH = Path("data/copper/slab/vacuum_layer")
FREE_LAYER_PATH = Path("data/copper/slab/free_layer")
SLAB_WIDTH_PATH = Path("data/copper/na/width")


def _prepare_vacuum_layer_convergence(atom: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for n_vacuum_layer in range(1, 9):
        slab_copper = get_surface(
            atom,
            (1, 1, 1),
            n_fixed_layer=5,
            n_vacuum_layer=n_vacuum_layer,
        )

        config = CastepConfig(
            data_path / f"slab_{n_vacuum_layer}_vacuum_layer",
            "slab",
        )
        params = SlabOptimizationParams(n_k_points=10, xc_functional="WC")
        calculator = get_slab_calculator(slab_copper, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


def _prepare_free_layer_convergence(atom: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)
    for n_free_layer in range(1, 6):
        n_fixed_layer = n_free_layer + 1
        slab_copper = get_surface(
            atom,
            (1, 1, 1),
            n_fixed_layer=n_fixed_layer,
            n_free_layer=n_free_layer,
            n_vacuum_layer=6,
        )

        config = CastepConfig(
            data_path / f"slab_{n_free_layer}_free_layer",
            "slab",
        )
        params = SlabOptimizationParams(n_k_points=12, xc_functional="WC")
        calculator = get_slab_calculator(
            slab_copper,
            params,
            config,
            task="GeometryOptimization",
        )
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


def _prepare_slab_width_convergence(
    atom: Atoms,
    data_path: Path,
) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    slab = get_surface(
        atom,
        (1, 1, 1),
        n_fixed_layer=4,
        n_free_layer=4,
        n_vacuum_layer=5,
    )

    slab_position = get_top_position(slab)
    rel_position = [0, 0, 1.0]
    na_atom = prepare_na_atom(slab_position, rel_position)

    for width in range(1, 9):
        adsorbate = prepare_adsorbate(na_atom, slab, width)
        config = CastepConfig(
            data_path / f"width_{width}",
            "slab_width",
        )
        params = SlabOptimizationParams(n_k_points=10, xc_functional="WC")
        calculator = get_slab_calculator(adsorbate, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


if __name__ == "__main__":
    bulk_config = CastepConfig(Path("data/copper/bulk/k_points_WC/bulk_10"), "bulk")
    bulk_copper = load_calculator_atoms(bulk_config)

    _prepare_free_layer_convergence(bulk_copper, FREE_LAYER_PATH)
    _prepare_vacuum_layer_convergence(bulk_copper, VACUUM_LAYER_PATH)
    _prepare_slab_width_convergence(bulk_copper, SLAB_WIDTH_PATH)

    slab_copper = get_surface(
        bulk_copper,
        (1, 1, 1),
        n_fixed_layer=5,
        n_vacuum_layer=5,
    )

    fig, _ = plot_atoms(slab_copper, radii=0.3, rotation=(-90, +30, 0))
    SAVE_DIR = Path(__file__).parent / "figures"
    plot_name = "initial_arrangement.png"
    save_fig(fig, SAVE_DIR / plot_name)
