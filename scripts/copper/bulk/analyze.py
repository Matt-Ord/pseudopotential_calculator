from pathlib import PosixPath

from pseudopotential_calculator.calculations.bulk import (
    plot_cell_length_against_cutoff_energy,
    plot_cell_length_against_n_k_points,
    plot_energy_against_cutoff_energy,
    plot_energy_against_n_k_points,
)
from pseudopotential_calculator.castep import load_all_calculators
from pseudopotential_calculator.scripting import maybe_copy_files_from_hpc

SAVE_DIR = PosixPath("scripts/copper/bulk")
K_POINTS_PATH = PosixPath("data/copper/bulk/k_points")
ENERGY_CUTOFF_PATH = PosixPath("data/copper/bulk/cutoff_energy")

if __name__ == "__main__":
    data_path = K_POINTS_PATH
    maybe_copy_files_from_hpc(data_path, K_POINTS_PATH)

    calculators = load_all_calculators(data_path)
    fig, ax, _ = plot_cell_length_against_n_k_points(calculators)

    # Add a horizontal line at y=2.53
    ax.axhline(y=2.53, color="red", linestyle="--", linewidth=1)  # type: ignore  # noqa: PGH003
    fig.show()
    fig.savefig(SAVE_DIR / "cell_length_against_n_k_points.png")  # type: ignore  # noqa: PGH003

    fig, _, _ = plot_energy_against_n_k_points(calculators)
    fig.tight_layout()
    fig.show()

    fig.savefig(SAVE_DIR / "energy_against_n_k_points")  # type: ignore

    data_path = ENERGY_CUTOFF_PATH
    maybe_copy_files_from_hpc(data_path, ENERGY_CUTOFF_PATH)

    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_cutoff_energy(calculators)
    fig.tight_layout()
    fig.show()

    fig.savefig(SAVE_DIR / "cell_length_against_cutoff_energy")  # type: ignore

    fig, ax, _ = plot_cell_length_against_cutoff_energy(calculators)
    # Add a horizontal line at y=2.53
    ax.axhline(y=2.53, color="red", linestyle="--", linewidth=1)  # type: ignore
    fig.tight_layout()
    fig.show()

    fig.savefig(SAVE_DIR / "energy_against_cutoff_energy")  # type: ignore

    input()
