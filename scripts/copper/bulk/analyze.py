from __future__ import annotations

from pathlib import Path, PosixPath

from pseudopotential_calculator.calculations.bulk import (
    plot_cell_length_against_cutoff_energy,
    plot_cell_length_against_n_k_points,
    plot_energy_against_cutoff_energy,
    plot_energy_against_n_k_points,
    plot_theoretical_cell_length,
)
from pseudopotential_calculator.castep import (
    CastepConfig,
    load_all_calculators,
    load_calculator_atoms,
)
from pseudopotential_calculator.scripting import maybe_copy_files_from_hpc
from pseudopotential_calculator.util import (
    get_figure,
    plot_atoms,
    save_fig,
    show_legend,
)

SAVE_DIR = Path(__file__).parent / "figures"

ENERGY_CUTOFF_PATH = PosixPath("data/copper/bulk/cutoff_energy")
K_POINTS_PATH_PBE = PosixPath("data/copper/bulk/k_points_PBE")
K_POINTS_PATH_WC = PosixPath("data/copper/bulk/k_points_WC")
K_POINTS_PATH_SP = PosixPath("data/copper/bulk/k_points_SP")


def _analyze_cell_length_convergence_with_n_k_points(paths: list[PosixPath]) -> None:
    fig, ax = get_figure(None)
    plot_theoretical_cell_length(ax, length=2.53)
    for data_path in paths:
        calculators = load_all_calculators(data_path)
        _, _, line = plot_cell_length_against_n_k_points(calculators, ax=ax)
        line.set_label(data_path.name)

    show_legend(ax)

    plot_name = "cell_length_against_n_k_points.png"
    save_fig(fig, SAVE_DIR / plot_name)


def _analyze_energy_convergence_with_n_k_points(paths: list[PosixPath]) -> None:
    fig, ax = get_figure(None)
    for data_path in paths:
        calculators = load_all_calculators(data_path)
        _, _, line = plot_energy_against_n_k_points(calculators, ax=ax)
        line.set_label(data_path.name)

    show_legend(ax)
    fig.tight_layout()
    plot_name = "energy_against_n_k_points.png"
    save_fig(fig, SAVE_DIR / plot_name)


def _analyze_cell_length_convergence_with_cutoff_energy(data_path: Path) -> None:
    calculators = load_all_calculators(data_path)

    fig, ax, _ = plot_cell_length_against_cutoff_energy(calculators)
    plot_theoretical_cell_length(ax, length=2.53)
    fig.tight_layout()
    show_legend(ax)

    plot_name = "cell_length_against_cutoff_energy.png"
    save_fig(fig, SAVE_DIR / plot_name)


def _analyze_energy_convergence_with_cutoff_energy(data_path: Path) -> None:
    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_cutoff_energy(calculators)
    fig.tight_layout()

    plot_name = "energy_against_cutoff_energy.png"
    save_fig(fig, SAVE_DIR / plot_name)


def _visualize_final_atoms(config: CastepConfig) -> None:
    atoms = load_calculator_atoms(config)

    fig, _ = plot_atoms(atoms, radii=0.3, rotation=(10, 0, 0))
    plot_name = "final_arrangement.png"
    save_fig(fig, SAVE_DIR / plot_name)


if __name__ == "__main__":
    for data_path in [
        K_POINTS_PATH_PBE,
        K_POINTS_PATH_WC,
        K_POINTS_PATH_SP,
        ENERGY_CUTOFF_PATH,
    ]:
        maybe_copy_files_from_hpc(data_path, data_path)

    _analyze_energy_convergence_with_cutoff_energy(ENERGY_CUTOFF_PATH)
    _analyze_cell_length_convergence_with_cutoff_energy(ENERGY_CUTOFF_PATH)

    k_points_paths = [
        K_POINTS_PATH_PBE,
        K_POINTS_PATH_WC,
        K_POINTS_PATH_SP,
    ]
    _analyze_energy_convergence_with_n_k_points(k_points_paths)
    _analyze_cell_length_convergence_with_n_k_points(k_points_paths)

    final_structure_path = Path("data/copper/bulk/k_points_WC/bulk_10")
    final_structure_config = CastepConfig(final_structure_path, "bulk")
    _visualize_final_atoms(final_structure_config)
