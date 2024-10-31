from __future__ import annotations

from pathlib import Path, PosixPath
from typing import TYPE_CHECKING

import numpy as np
from ase import Atom

from pseudopotential_calculator.atoms import Adsorbate, add_vectors
from pseudopotential_calculator.calculations.adsorbate import (
    append_atom_to_atoms,
    get_top_position,
    prepare_adsorbate_positions,
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


adsorbate_config = Adsorbate(name="sodium", symbol="Na", height=2.0, cut_off_energy=800)

ENERGY_MAP_PATH = Path(f"data/copper/{adsorbate_config.name}/energy_map")
VACUUM_HEIGHT_PATH = Path(f"data/copper/{adsorbate_config.name}/vacuum_height")
SLAB_WIDTH_PATH = Path(f"data/copper/{adsorbate_config.name}/width")


def _prepare_vacuum_height_convergence(atoms: Atoms, data_path: Path) -> None:
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for vacuum_height in range(1, 10):
        slab = get_surface(
            atoms,
            (1, 1, 1),
            n_fixed_layer=4,
            n_free_layer=3,
            vacuum_height=vacuum_height + adsorbate_config.height,
        )

        slab_position = get_top_position(slab)
        adsorbate_position = add_vectors(
            slab_position,
            (0, 0, adsorbate_config.height),
        )
        adsorbate_atom = Atom(adsorbate_config.symbol, adsorbate_position)

        slab_with_adsorbate = append_atom_to_atoms(adsorbate_atom, slab)

        config = CastepConfig(
            data_path / f"{vacuum_height}",
            f"{adsorbate_config.symbol}_{vacuum_height}_vacuum",
        )
        params = CastepParams(
            n_k_points=(12, 12, 1),
            xc_functional="WC",
            cut_off_energy=adsorbate_config.cut_off_energy,
        )
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

    slab = get_surface(
        atoms,
        (1, 1, 1),
        n_fixed_layer=4,
        n_free_layer=3,
        vacuum_height=10 + adsorbate_config.height,
    )

    slab_position = get_top_position(slab)
    adsorbate_position = add_vectors(
        slab_position,
        (0, 0, adsorbate_config.height),
    )
    adsorbate_atom = Atom(adsorbate_config.symbol, adsorbate_position)

    for width in range(1, 8):
        slab_with_adsorbate = prepare_slab_with_adsorbate(adsorbate_atom, slab, width)
        config = CastepConfig(
            data_path / f"{width}",
            f"{adsorbate_config.symbol}_{width}_width",
        )
        params = CastepParams(
            n_k_points=(int(np.ceil(12 / width)), int(np.ceil(12 / width)), 1),
            xc_functional="WC",
            max_scf_cycles=100,
            cut_off_energy=adsorbate_config.cut_off_energy,
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
    adsorbate_positions = prepare_adsorbate_positions(
        slab,
        heights_above_surface,
        scale_factor=scale_factor,
    )
    for position in adsorbate_positions:
        adsorbate_atom = Atom(adsorbate_config.symbol, position)
        adsorbate = prepare_slab_with_adsorbate(adsorbate_atom, slab, width)

        config = CastepConfig(
            data_path / f"{round(position[0], 1)}x{round(position[1], 1)}x"
            f"{round(position[2], 1)}",
            f"{adsorbate_config.symbol}_energy_map",
        )
        params = CastepParams(
            n_k_points=(12, 12, 1),
            xc_functional="WC",
            cut_off_energy=adsorbate_config.cut_off_energy,
        )
        calculator = get_calculator(adsorbate, params, config)
        prepare_calculator_with_submit_script(calculator)
        calculators.append(calculator)

    prepare_submit_all_script(calculators, data_path)
    maybe_copy_files_to_hpc(data_path, PosixPath(data_path.as_posix()))


if __name__ == "__main__":
    bulk_config = CastepConfig(Path("data/copper/bulk/k_points_WC/bulk_10"), "bulk")
    bulk_copper = load_calculator_atoms(bulk_config)

    _prepare_vacuum_height_convergence(bulk_copper, VACUUM_HEIGHT_PATH)
    _prepare_slab_width_convergence(bulk_copper, SLAB_WIDTH_PATH)

    slab_config = CastepConfig(
        Path("data/copper/slab/free_layer/slab_4_free_layer"),
        "slab_free_layer",
    )
    slab_copper = load_calculator_atoms(slab_config)

    heights_above_surface = [2.0, 5.0]
    _prepare_energy_map(
        slab_copper,
        ENERGY_MAP_PATH,
        heights_above_surface,
        width=2,
        scale_factor=6,
    )
