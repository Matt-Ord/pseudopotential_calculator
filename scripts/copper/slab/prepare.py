from pathlib import Path, PosixPath

from ase import Atoms
from ase.build import (
    bulk,
    surface,  # type: ignore  # noqa: PGH003
)
from ase.io import read  # type: ignore  # noqa: PGH003

from pseudopotential_calculator.calculations.slab import (
    SlabOptimizationParams,
    get_slab_optimization_calculator,
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

    for n_vaccum_layer in range(3, 21):
        slab_copper = surface(atom, (1, 1, 0), 5, vacuum=n_vaccum_layer)  # type: ignore
        config = CastepConfig(
            data_path / f"slab_{n_vaccum_layer}_vaccum_layer",
            "slab",
        )
        params = SlabOptimizationParams(n_k_points=10)
        calculator = get_slab_optimization_calculator(slab_copper, params, config)  # type: ignore
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


VACUUM_LAYER_PATH = Path("data/copper/slab/vaccum_layer")  # type: ignore

if __name__ == "__main__":
    bulk_copper = read("Cu_test.poscar", format="vasp")

    bulk_copper = bulk("Cu", "fcc", 3.57743)
    # TODO: PUTS IN Wrong directory

    # write the model file to be viewed
    _prepare_vaccum_layer_convergence(bulk_copper, VACUUM_LAYER_PATH)  # type: ignore
