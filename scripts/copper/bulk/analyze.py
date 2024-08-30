from pathlib import Path, PosixPath

from pseudopotential_calculator.calculations.bulk import (
    plot_cell_length_against_cutoff_energy,
    plot_cell_length_against_n_k_points,
    plot_energy_against_cutoff_energy,
    plot_energy_against_n_k_points,
)
from pseudopotential_calculator.castep import load_all_calculators
from pseudopotential_calculator.scripting import maybe_copy_files_from_hpc
from pseudopotential_calculator.util import get_figure

SAVE_DIR = Path(__file__).parent / "figures"

ENERGY_CUTOFF_PATH = PosixPath("data/copper/bulk/cutoff_energy")
K_POINTS_PATH_PBE = PosixPath("data/copper/bulk/k_points_PBE")
K_POINTS_PATH_WC = PosixPath("data/copper/bulk/k_points_WC")
K_POINTS_PATH_SP = PosixPath("data/copper/bulk/k_points_SP")

if __name__ == "__main__":
    paths = [
        ENERGY_CUTOFF_PATH,
        K_POINTS_PATH_PBE,
        K_POINTS_PATH_WC,
        K_POINTS_PATH_SP,
    ]

    fig1, ax1 = get_figure(None)
    ax1.axhline(y=2.53, linestyle="--")  # type: ignore  # noqa: PGH003

    for data_path in paths:
        maybe_copy_files_from_hpc(data_path, data_path)
        calculators = load_all_calculators(data_path)
        _, _, line = plot_cell_length_against_n_k_points(calculators, ax=ax1)
        line.set_label(data_path.name)

    ax1.legend()
    fig1.show()
    fig1.savefig(SAVE_DIR / "cell_length_against_n_k_points.Png")  # type: ignore  # noqa: PGH003
    data_path = K_POINTS_PATH_PBE

    fig2, ax2 = get_figure(None)
    for data_path in paths:
        maybe_copy_files_from_hpc(data_path, data_path)
        calculators = load_all_calculators(data_path)
        _, _, line = plot_energy_against_n_k_points(calculators, ax=ax2)
        line.set_label(data_path.name)

    ax2.legend()
    fig2.show()
    fig2.savefig(SAVE_DIR / "energy_against_n_k_points.Png")  # type: ignore  # noqa: PGH003

    fig3, ax3 = get_figure(None)
    data_path = ENERGY_CUTOFF_PATH
    maybe_copy_files_from_hpc(data_path, ENERGY_CUTOFF_PATH)

    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_cutoff_energy(calculators)

    fig.tight_layout()
    fig.show()
    fig.savefig(SAVE_DIR / "cell_length_against_cutoff_energy")  # type: ignore

    fig, ax, _ = plot_cell_length_against_cutoff_energy(calculators)
    ax.axhline(y=2.53, linestyle="--")  # type: ignore
    fig.tight_layout()
    fig.show()

    fig.savefig(SAVE_DIR / "energy_against_cutoff_energy")  # type: ignore
