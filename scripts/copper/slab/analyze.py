from pathlib import Path, PosixPath

from pseudopotential_calculator.calculations.slab import (
    plot_displacement_against_n_free_layer,
    plot_energy_against_n_vacuum_layer,
)
from pseudopotential_calculator.castep import (
    load_all_calculators,
)
from pseudopotential_calculator.scripting import maybe_copy_files_from_hpc
from pseudopotential_calculator.util import get_figure, save_fig, show_legend

VACUUM_LAYER_PATH = PosixPath("data/copper/slab/vacuum_layer")
FREE_LAYER_PATH = PosixPath("data/copper/slab/free_layer")

SAVE_DIR = Path(__file__).parent / "figures"


def _analyze_convergence_with_n_vacuum_layer() -> None:
    data_path = VACUUM_LAYER_PATH
    maybe_copy_files_from_hpc(data_path, data_path)
    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_n_vacuum_layer(calculators)
    fig.tight_layout()
    plot_name = "energy_against_n_vacuum_layer.png"
    save_fig(fig, SAVE_DIR / plot_name)


def _analyze_convergence_with_n_free_layer() -> None:
    data_path = FREE_LAYER_PATH
    maybe_copy_files_from_hpc(data_path, data_path)
    calculators = load_all_calculators(data_path)
    fig, ax = get_figure(None)
    for layer in range(6):
        _, _, line = plot_displacement_against_n_free_layer(
            calculators,
            -(layer + 1),
            ax=ax,
        )
        line.set_label(f"layer {layer}")
    show_legend(ax)
    plot_name = "displacement_against_n_free_layer.png"
    save_fig(fig, SAVE_DIR / plot_name)


if __name__ == "__main__":
    _analyze_convergence_with_n_free_layer()
    _analyze_convergence_with_n_vacuum_layer()
