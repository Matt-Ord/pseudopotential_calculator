from pathlib import Path, PosixPath

import numpy as np
from ase import Atoms
from ase.build import (
    bulk,  # type: ignore bad library
)

from pseudopotential_calculator.calculations.slab import (
    SlabOptimizationParams,
    get_slab_vaccum_calculator,
    get_surface,
)
from pseudopotential_calculator.castep import Castep, CastepConfig
from pseudopotential_calculator.hpc import (
    copy_files_to_hpc,
    prepare_calculator_with_submit_script,
    prepare_submit_all_script,
)
from pseudopotential_calculator.util import prepare_clean_directory


def _prepare_vaccum_layer_convergence(atom: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)
    height_per_vaccum_layer = 1.788715
    for n_vaccum_layer in range(3, 21):
        slab_copper = get_surface(
            atom,
            (1, 1, 1),
            5,
            vacuum=n_vaccum_layer * height_per_vaccum_layer,
        )

        zmin = np.min(slab_copper.positions[:, 2])  # type: ignore

        # Adjust the positions to move the bottom layer to z = 0
        slab_copper.positions[:, 2] -= zmin  # type: ignore

        config = CastepConfig(
            data_path / f"slab_{n_vaccum_layer}_vaccum_layer",
            "slab",
        )
        params = SlabOptimizationParams(n_k_points=10)
        calculator = get_slab_vaccum_calculator(slab_copper, params, config)  # type: ignore
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


VACUUM_LAYER_PATH = Path("data/copper/slab/vaccum_layer")

if __name__ == "__main__":
    bulk_copper = bulk("Cu", "fcc", 3.57743)
    # TODO: PUTS IN Wrong directory
    _prepare_vaccum_layer_convergence(bulk_copper, VACUUM_LAYER_PATH)  # type: ignore
