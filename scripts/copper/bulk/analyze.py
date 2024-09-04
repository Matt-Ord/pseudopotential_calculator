from __future__ import annotations

from pathlib import Path, PosixPath
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from ase.build import (
    bulk,  # type: ignore bad library
)
from ase.visualize.plot import plot_atoms  # type: ignore bad lib

from pseudopotential_calculator.calculations.bulk import (
    plot_cell_length_against_cutoff_energy,
    plot_cell_length_against_n_k_points,
    plot_energy_against_cutoff_energy,
    plot_energy_against_n_k_points,
    plot_theoretical_cell_length,
)
from pseudopotential_calculator.castep import (
    CastepConfig,
    get_calculator_atom,
    load_all_calculators,
    load_calculator,
)
from pseudopotential_calculator.scripting import maybe_copy_files_from_hpc
from pseudopotential_calculator.util import get_figure, save_fig, show_legend

if TYPE_CHECKING:
    from ase import Atoms

SAVE_DIR = Path(__file__).parent / "figures"

ENERGY_CUTOFF_PATH = PosixPath("data/copper/bulk/cutoff_energy")
K_POINTS_PATH_PBE = PosixPath("data/copper/bulk/k_points_PBE")
K_POINTS_PATH_WC = PosixPath("data/copper/bulk/k_points_WC")
K_POINTS_PATH_SP = PosixPath("data/copper/bulk/k_points_SP")


def _download_all() -> None:
    download_paths = [
        ENERGY_CUTOFF_PATH,
        K_POINTS_PATH_PBE,
        K_POINTS_PATH_WC,
        K_POINTS_PATH_SP,
    ]

    for data_path in download_paths:
        maybe_copy_files_from_hpc(data_path, data_path)


def _analyze_cell_length_convergence_with_n_k_points(paths: list[PosixPath]) -> None:
    fig, ax = get_figure(None)
    plot_theoretical_cell_length(ax, length=2.53)
    for data_path in paths:
        calculators = load_all_calculators(data_path)
        _, _, line = plot_cell_length_against_n_k_points(calculators, ax=ax)
        line.set_label(data_path.name)

    show_legend(ax)
    fig.show()
    plot_name = "cell_length_against_n_k_points.Png"
    save_fig(fig, SAVE_DIR / plot_name)


def _analyze_energy_convergence_with_n_k_points(paths: list[PosixPath]) -> None:
    fig, ax = get_figure(None)
    for data_path in paths:
        calculators = load_all_calculators(data_path)
        _, _, line = plot_energy_against_n_k_points(calculators, ax=ax)
        line.set_label(data_path.name)
    show_legend(ax)
    fig.show()
    plot_name = "energy_against_n_k_points.Png"
    save_fig(fig, SAVE_DIR / plot_name)


def _analyze_cell_length_convergence_with_cutoff_energy() -> None:
    data_path = ENERGY_CUTOFF_PATH

    calculators = load_all_calculators(data_path)

    fig, ax, _ = plot_cell_length_against_cutoff_energy(calculators)
    plot_theoretical_cell_length(ax, length=2.53)
    fig.tight_layout()
    show_legend(ax)
    fig.show()
    plot_name = "cell_length_against_cutoff_energy"
    save_fig(fig, SAVE_DIR / plot_name)


def _analyze_energy_convergence_with_cutoff_energy() -> None:
    data_path = ENERGY_CUTOFF_PATH
    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_cutoff_energy(calculators)
    fig.tight_layout()
    fig.show()
    plot_name = "energy_against_cutoff_energy"
    save_fig(fig, SAVE_DIR / plot_name)


def _visualize_initial_cell(atom: Atoms | None) -> None:
    fig, ax = plt.subplots()  # type: ignore bad library
    plot_atoms(atom, ax, radii=0.3, rotation=("10x,0y,0z"))
    plot_name = "initial_arrangement"
    save_fig(fig, SAVE_DIR / plot_name)


def _visualize_final_cell(atom: Atoms | None) -> None:
    if atom is None:
        ValueError("No atom")
    fig, ax = plt.subplots()  # type: ignore bad library
    plot_atoms(atom, ax, radii=0.3, rotation=("10x,0y,0z"))
    plot_name = "final_arrangement"
    save_fig(fig, SAVE_DIR / plot_name)


if __name__ == "__main__":
    paths = [
        K_POINTS_PATH_PBE,
        K_POINTS_PATH_WC,
        K_POINTS_PATH_SP,
    ]
    initial_structure = bulk("Cu", "fcc", 3.8)

    final_structure_path = Path("data/copper/bulk/k_points_WC/bulk_10")
    final_structure_config = CastepConfig(final_structure_path, "bulk")
    final_structure = get_calculator_atom(load_calculator(final_structure_config))
    _download_all()
    _visualize_initial_cell(initial_structure)
    _visualize_final_cell(final_structure)
    _analyze_energy_convergence_with_cutoff_energy()
    _analyze_cell_length_convergence_with_cutoff_energy()
    _analyze_energy_convergence_with_n_k_points(paths)
    _analyze_cell_length_convergence_with_n_k_points(paths)
