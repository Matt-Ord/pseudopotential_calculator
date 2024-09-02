from __future__ import annotations

from pathlib import Path, PosixPath
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from ase.build import (
    bulk,  # type: ignore bad library
)
from ase.visualize.plot import plot_atoms  # type: ignore

from pseudopotential_calculator.calculations.bulk import (
    plot_cell_length_against_cutoff_energy,
    plot_cell_length_against_n_k_points,
    plot_energy_against_cutoff_energy,
    plot_energy_against_n_k_points,
    plot_theoretical_cell_length,
)
from pseudopotential_calculator.castep import load_all_calculators
from pseudopotential_calculator.scripting import maybe_copy_files_from_hpc
from pseudopotential_calculator.util import get_figure, save_fig, show_legend

if TYPE_CHECKING:
    from ase import Atoms

SAVE_DIR = Path(__file__).parent / "figures"

ENERGY_CUTOFF_PATH = PosixPath("data/copper/bulk/cutoff_energy")
K_POINTS_PATH_PBE = PosixPath("data/copper/bulk/k_points_PBE")
K_POINTS_PATH_WC = PosixPath("data/copper/bulk/k_points_WC")
K_POINTS_PATH_SP = PosixPath("data/copper/bulk/k_points_SP")


def _plot_cell_length_against_n_k_points(paths: list[PosixPath]) -> None:
    fig, ax = get_figure(None)
    plot_theoretical_cell_length(ax)
    for data_path in paths:
        maybe_copy_files_from_hpc(data_path, data_path)
        calculators = load_all_calculators(data_path)
        _, _, line = plot_cell_length_against_n_k_points(calculators, ax=ax)
        line.set_label(data_path.name)

    show_legend(ax)
    fig.show()
    plot_name = "cell_length_against_n_k_points.Png"
    save_fig(fig, SAVE_DIR / plot_name)


def _plot_energy_against_n_k_points(paths: list[PosixPath]) -> None:
    fig, ax = get_figure(None)
    for data_path in paths:
        maybe_copy_files_from_hpc(data_path, data_path)
        calculators = load_all_calculators(data_path)
        _, _, line = plot_energy_against_n_k_points(calculators, ax=ax)
        line.set_label(data_path.name)
    show_legend(ax)
    fig.show()
    plot_name = "energy_against_n_k_points.Png"
    save_fig(fig, SAVE_DIR / plot_name)


def _plot_cell_length_against_cutoff_energy() -> None:
    data_path = ENERGY_CUTOFF_PATH
    maybe_copy_files_from_hpc(data_path, ENERGY_CUTOFF_PATH)
    calculators = load_all_calculators(data_path)

    fig, ax, _ = plot_cell_length_against_cutoff_energy(calculators)
    plot_theoretical_cell_length(ax)
    fig.tight_layout()
    show_legend(ax)
    fig.show()
    plot_name = "cell_length_against_cutoff_energy"
    save_fig(fig, SAVE_DIR / plot_name)


def _plot_energy_against_cutoff_energy() -> None:
    data_path = ENERGY_CUTOFF_PATH
    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_cutoff_energy(calculators)
    fig.tight_layout()
    fig.show()
    plot_name = "energy_against_cutoff_energy"
    save_fig(fig, SAVE_DIR / plot_name)


def _visualize_initial_cell(atom: Atoms) -> None:
    fig, ax = plt.subplots()  # type: ignore bad library
    plot_atoms(atom, ax, radii=0.3, rotation=("45x,45y,0z"))
    plot_name = "initial arrangement"
    save_fig(fig, SAVE_DIR / plot_name)


def _visualize_final_cell(atom: Atoms) -> None:
    fig, ax = plt.subplots()  # type: ignore bad library
    plot_atoms(atom, ax, radii=0.3, rotation=("0x,0y,0z"))
    plot_name = "final arrangement"
    save_fig(fig, SAVE_DIR / plot_name)


if __name__ == "__main__":
    paths = [
        ENERGY_CUTOFF_PATH,
        K_POINTS_PATH_PBE,
        K_POINTS_PATH_WC,
        K_POINTS_PATH_SP,
    ]

    initial_structure = bulk("Cu", "fcc", 3.8)
    # TODO load final structre
    final_structure = bulk("Cu", "fcc", 3.8)
    _visualize_initial_cell(initial_structure)
    _visualize_final_cell(final_structure)
    _plot_energy_against_cutoff_energy()
    _plot_cell_length_against_cutoff_energy()
    _plot_energy_against_n_k_points(paths)
    _plot_cell_length_against_n_k_points(paths)
