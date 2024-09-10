from pathlib import Path, PosixPath

from ase import Atoms

from pseudopotential_calculator.calculations.adsorbate import (
    get_basis_vectors,
    get_top_layer_z_vector,
    prepare_adsorbate,
    prepare_na_atom,
)
from pseudopotential_calculator.calculations.slab import (
    SlabOptimizationParams,
    get_slab_calculator,
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

ENERGY_MAP_PATH = Path("data/copper/energy_map")


def prepare_energy_map(
    slab: Atoms,
    width: int,
    data_path: Path,
    n_division: int,
    n_layer: int,
    resolution_in_z: float,
) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    basis_vectors = get_basis_vectors(slab)
    z = get_top_layer_z_vector(slab)
    x = basis_vectors[0]
    y = basis_vectors[1]

    for i in range(n_division + 1):
        for j in range(i + 1):
            slab_position = [i * x + j * y + z for x, y, z in zip(x, y, z)]

            for k in range(n_layer + 1):
                rel_position = [0, 0, k * resolution_in_z]
                na_atom = prepare_na_atom(slab_position, rel_position)
                adsorbate = prepare_adsorbate(na_atom, slab, width)

                config = CastepConfig(
                    data_path / f"{i}x{j}",
                    "energy_map",
                )
                params = SlabOptimizationParams(n_k_points=10, xc_functional="WC")
                calculator = get_slab_calculator(adsorbate, params, config)
                prepare_calculator_with_submit_script(calculator)
                calculators.append(calculator)

    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


if __name__ == "__main__":
    slab_config = CastepConfig(
        Path("data/copper/slab/free_layer/slab_3_free_layer"),
        "slab",
    )
    slab = load_calculator_atoms(slab_config)
    prepare_energy_map(slab, 3, ENERGY_MAP_PATH, 6, 1, 2)
