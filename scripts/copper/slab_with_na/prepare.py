from pathlib import Path, PosixPath

from ase import Atom, Atoms

from pseudopotential_calculator.calculations.slab import (
    SlabOptimizationParams,
    add_atom_onto_slab,
    get_slab_calculator,
    get_top_position,
    repeat_slab,
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

SLAB_SIZE_PATH = Path("data/copper/slab/slab_size")
FREE_LAYER_PATH = Path("data/copper/slab/free_layer")


def _prepare_na_atom(slab_copper: Atoms) -> Atom:
    na_rel_position = [0, 0, 1]
    cu_position = get_top_position(slab_copper)
    na_position = [x + y for x, y in zip(cu_position, na_rel_position)]
    return Atom("Na", position=na_position)


def _prepare_slab_size_convergence(
    slab_copper: Atoms,
    na_atom: Atom,
    data_path: Path,
) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for n in range(1, 9):
        repeated_slab = repeat_slab(slab_copper, (n, n, 1))
        add_atom_onto_slab(na_atom, repeated_slab)
        config = CastepConfig(
            data_path / f"slab_{n}_repeats",
            "slab",
        )
        params = SlabOptimizationParams(n_k_points=10, xc_functional="WC")
        calculator = get_slab_calculator(repeated_slab, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


if __name__ == "__main__":
    # TODO change to correct path to read file when data is pulled from hpc
    slab_config = CastepConfig(
        Path("data/copper/slab/vaccum_layer/slab_5_vaccum_layer"),
        "slab",
    )
    slab_copper = load_calculator_atoms(slab_config)

    na_atom = _prepare_na_atom(slab_copper)
    _prepare_slab_size_convergence(slab_copper, na_atom, SLAB_SIZE_PATH)
