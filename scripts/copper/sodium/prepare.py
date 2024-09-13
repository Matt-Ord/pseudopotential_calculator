from __future__ import annotations

from pathlib import Path, PosixPath
from typing import TYPE_CHECKING

import numpy as np

from pseudopotential_calculator.calculations.adsorbate import (
    get_adsorbate_position,
    get_scaled_basis_vectors,
    get_slab_with_adsorbate,
    get_top_layer_z_vector,
    get_top_position,
    prepare_adsorbate,
    prepare_slab_with_adsorbate,
)
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
    prepare_clean_directory,
)

if TYPE_CHECKING:
    from ase import Atoms

ENERGY_MAP_PATH = Path("data/copper/sodium/energy_map")
VACUUM_LAYER_PATH = Path("data/copper/sodium/vacuum_layer")
SLAB_WIDTH_PATH = Path("data/copper/sodium/width")


def _prepare_vacuum_layer_convergence(atoms: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    adsorbate_height = 1.0

    for n_vacuum_layer in range(1, 9):
        slab = get_surface(
            atoms,
            (1, 1, 1),
            n_fixed_layer=4,
            n_free_layer=3,
            n_vacuum_layer=n_vacuum_layer,
            adsorbate_height=adsorbate_height,
        )

        slab_position = get_top_position(slab)
        adsorbate_position = get_adsorbate_position(
            slab_position,
            adsorbate_height,
        )
        na_atom = prepare_adsorbate("Na", adsorbate_position)

        slab_with_adsorbate = get_slab_with_adsorbate(na_atom, slab)

        config = CastepConfig(
            data_path / f"slab_{n_vacuum_layer}_vacuum_layer",
            "slab",
        )
        params = CastepParams(n_k_points=(12, 12, 1), xc_functional="WC")
        calculator = get_calculator(slab_with_adsorbate, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


def _prepare_slab_width_convergence(
    atoms: Atoms,
    data_path: Path,
) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)
    adsorbate_height = 1.0

    slab = get_surface(
        atoms,
        (1, 1, 1),
        n_fixed_layer=4,
        n_free_layer=3,
        n_vacuum_layer=5,
        adsorbate_height=adsorbate_height,
    )

    slab_position = get_top_position(slab)
    adsorbate_position = get_adsorbate_position(
        slab_position,
        adsorbate_height,
    )
    na_atom = prepare_adsorbate("Na", adsorbate_position)

    for width in range(1, 9):
        slab_with_adsorbate = prepare_slab_with_adsorbate(na_atom, slab, width)
        config = CastepConfig(
            data_path / f"width_{width}",
            "slab_width",
        )
        params = CastepParams(
            n_k_points=(int(np.ceil(12 / width)), int(np.ceil(12 / width)), 1),
            xc_functional="WC",
            max_scf_cycles=100,
        )
        calculator = get_calculator(slab_with_adsorbate, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)
    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


def _prepare_energy_map(
    slab: Atoms,
    data_path: Path,
    heights_above_surface: list[float],
    *,
    width: int,
    scale_factor: int,
) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    basis_vectors = get_scaled_basis_vectors(slab, scale_factor=scale_factor)
    z = get_top_layer_z_vector(slab)
    x = basis_vectors[0]
    y = basis_vectors[1]

    i_vals, j_vals, k_vals = np.meshgrid(
        np.arange(scale_factor + 1),
        np.arange(scale_factor + 1),
        heights_above_surface,
        indexing="ij",
    )

    i_flat = i_vals.ravel()
    j_flat = j_vals.ravel()
    k_flat = k_vals.ravel()

    for idx in range(len(i_flat)):
        i, j, k = i_flat[idx], j_flat[idx], k_flat[idx]

        if i >= j:
            slab_position = (
                i * x[0] + j * y[0] + z[0],
                i * x[1] + j * y[1] + z[1],
                i * x[2] + j * y[2] + z[2],
            )

            adsorbate_position = get_adsorbate_position(slab_position, k)
            na_atom = prepare_adsorbate("Na", adsorbate_position)
            adsorbate = prepare_slab_with_adsorbate(na_atom, slab, width)

            config = CastepConfig(
                data_path / f"{i}x{j}x{k}",
                "energy_map",
            )
            params = CastepParams(n_k_points=(12, 12, 1), xc_functional="WC")
            calculator = get_calculator(adsorbate, params, config)
            prepare_calculator_with_submit_script(calculator)
            calculators.append(calculator)

    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


if __name__ == "__main__":
    bulk_config = CastepConfig(Path("data/copper/bulk/k_points_WC/bulk_10"), "bulk")
    bulk_copper = load_calculator_atoms(bulk_config)

    _prepare_vacuum_layer_convergence(bulk_copper, VACUUM_LAYER_PATH)
    _prepare_slab_width_convergence(bulk_copper, SLAB_WIDTH_PATH)

    slab_config = CastepConfig(
        Path("data/copper/slab/free_layer/slab_3_free_layer"),
        "slab",
    )
    slab_copper = load_calculator_atoms(slab_config)
    heights_above_surface = [1.0, 5.0]
    _prepare_energy_map(
        slab_copper,
        ENERGY_MAP_PATH,
        heights_above_surface,
        width=3,
        scale_factor=6,
    )
