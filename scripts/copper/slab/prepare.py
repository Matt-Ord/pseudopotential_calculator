from pathlib import Path, PosixPath
from typing import cast

from ase import Atoms

from pseudopotential_calculator.calculations.slab import (
    SlabOptimizationParams,
    get_slab_vaccum_calculator,
    get_surface,
)
from pseudopotential_calculator.castep import (
    Castep,
    CastepConfig,
    get_calculator_atom,
    load_calculator,
)
from pseudopotential_calculator.hpc import (
    prepare_calculator_with_submit_script,
    prepare_submit_all_script,
)
from pseudopotential_calculator.scripting import maybe_copy_files_to_hpc
from pseudopotential_calculator.util import prepare_clean_directory


def _prepare_vaccum_layer_convergence(atom: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for n_vaccum_layer in range(3, 21):
        slab_copper = get_surface(
            atom,
            (1, 1, 1),
            5,
            n_vaccum_layer,
        )

        config = CastepConfig(
            data_path / f"slab_{n_vaccum_layer}_vaccum_layer",
            "slab",
        )
        params = SlabOptimizationParams(n_k_points=10)
        calculator = get_slab_vaccum_calculator(slab_copper, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


VACUUM_LAYER_PATH = Path("data/copper/slab/vaccum_layer")

if __name__ == "__main__":
    data_path = Path("data/copper/bulk/k_points_WC/bulk_10")
    config = CastepConfig(data_path, "bulk")
    bulk_copper = cast(Atoms, get_calculator_atom(load_calculator(config)))
    _prepare_vaccum_layer_convergence(bulk_copper, VACUUM_LAYER_PATH)
