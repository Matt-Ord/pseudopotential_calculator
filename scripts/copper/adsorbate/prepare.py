from pathlib import Path, PosixPath

from ase import Atoms

from pseudopotential_calculator.calculations.adsorbate import (
    AdsorbateOptimizationParams,
    get_adsorbate_optimization_calculator,
    prepare_adsorbate,
    prepare_na_atom,
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
    prepare_clean_directory,
)

SLAB_WIDTH_PATH = Path("data/copper/na/width")
FREE_LAYER_PATH = Path("data/copper/slab/free_layer")


def _prepare_slab_size_convergence(
    slab_copper: Atoms,
    data_path: Path,
) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)
    na_atom = prepare_na_atom(slab_copper)
    for width in range(1, 9):
        adsorbate = prepare_adsorbate(na_atom, slab_copper, width)
        config = CastepConfig(
            data_path / f"width_{width}",
            "slab_width",
        )
        params = AdsorbateOptimizationParams(n_k_points=10, xc_functional="WC")
        calculator = get_adsorbate_optimization_calculator(adsorbate, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


if __name__ == "__main__":
    slab_config = CastepConfig(
        Path("data/copper/slab/free_layer/slab_3_free_layer"),
        "slab",
    )
    slab_copper = load_calculator_atoms(slab_config)

    _prepare_slab_size_convergence(slab_copper, SLAB_WIDTH_PATH)
