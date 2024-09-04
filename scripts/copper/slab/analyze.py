from pathlib import Path, PosixPath

import matplotlib.pyplot as plt
from ase import Atoms
from ase.visualize.plot import plot_atoms  # type: ignore bad lib

from pseudopotential_calculator.calculations.slab import (
    get_surface,
    plot_energy_against_n_vacuum_layer,
)
from pseudopotential_calculator.castep import (
    CastepConfig,
    get_calculator_atom,
    load_all_calculators,
    load_calculator,
)
from pseudopotential_calculator.scripting import maybe_copy_files_from_hpc
from pseudopotential_calculator.util import save_fig

VACUUM_LAYER_PATH = PosixPath("data/copper/slab/vaccum_layer")

SAVE_DIR = Path(__file__).parent / "figures"


def _visualize_initial_slab(atom: Atoms) -> None:
    fig, ax = plt.subplots()  # type: ignore bad library
    plot_atoms(atom, ax, radii=0.3, rotation=("10x,0y,0z"))
    plot_name = "initial_arrangement"
    save_fig(fig, SAVE_DIR / plot_name)


def _analyze_convergence_with_n_vaccum_layer() -> None:
    data_path = VACUUM_LAYER_PATH
    maybe_copy_files_from_hpc(data_path, data_path)
    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_n_vacuum_layer(calculators)
    fig.tight_layout()
    plot_name = "energy_against_n_vaccum_layer"
    save_fig(fig, SAVE_DIR / plot_name)


def _get_initial_slab(data_path: Path) -> Atoms:
    config = CastepConfig(data_path, "bulk")
    atom = get_calculator_atom(load_calculator(config))
    if atom is None:
        msg = "The 'atom' parameter cannot be None."
        raise ValueError(msg)
    return get_surface(
        atom,
        (1, 1, 1),
        5,
        n_vacuum_layer=0,
    )


if __name__ == "__main__":
    bulk_data_path = Path("data/copper/bulk/k_points_WC/bulk_10")
    initial_slab = _get_initial_slab(bulk_data_path)
    _visualize_initial_slab(initial_slab)
    _analyze_convergence_with_n_vaccum_layer()
