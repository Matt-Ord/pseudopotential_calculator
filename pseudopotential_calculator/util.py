from __future__ import annotations

from typing import TYPE_CHECKING, cast

import matplotlib.pyplot as plt
import numpy as np
from ase.visualize.plot import plot_atoms as plot_atoms_ase  # type: ignore bad lib
from matplotlib.axes import Axes
from matplotlib.figure import Figure

if TYPE_CHECKING:
    from pathlib import Path

    from ase import Atoms
    from matplotlib.lines import Line2D


def _rename_calculation_folder(directory: Path) -> None:
    if directory.exists():
        n = 1
        while (directory.parent / f"{directory.name}.old.{n}").exists():
            n += 1
        directory.rename(directory.parent / f"{directory.name}.old.{n}")


def prepare_clean_directory(directory: Path) -> None:
    _rename_calculation_folder(directory)
    directory.mkdir(parents=True, exist_ok=True)


def get_figure(ax: Axes | None) -> tuple[Figure, Axes]:
    """Get the figure of the given axis.

    If no figure exists, a new figure is created

    Parameters
    ----------
    ax : Axes | None

    Returns
    -------
    tuple[Figure, Axes]

    """
    if ax is None:
        return cast(tuple[Figure, Axes], plt.subplots())  # type: ignore plt.subplots Unknown type

    fig = ax.get_figure()
    if fig is None:
        fig = plt.figure()  # type: ignore plt.figure Unknown type
        ax.set_figure(fig)
    return fig, ax


def _get_data_label(name: str, unit: str | None) -> str:
    if unit is None:
        return name
    return f"{name} / {unit}"


def plot_data_comparison(
    x: tuple[str, np.ndarray[tuple[int], np.dtype[np.float64]], str | None],
    y: tuple[str, np.ndarray[tuple[int], np.dtype[np.float64]], str | None],
    *,
    ax: Axes | None = None,
) -> tuple[Figure, Axes, Line2D]:
    fig, ax = get_figure(ax)
    x_label, x_data, x_unit = x
    y_label, y_data, y_unit = y
    args = np.argsort(x_data)
    (line,) = ax.plot(x_data[args], y_data[args])  # type: ignore bad library

    ax.set_xlabel(_get_data_label(x_label, x_unit))  # type: ignore bad library
    ax.set_ylabel(_get_data_label(y_label, y_unit))  # type: ignore bad library
    ax.set_title(f"{y_label} vs {x_label}")  # type: ignore bad library

    return fig, ax, line


def plot_atoms(
    atoms: Atoms,
    *,
    ax: Axes | None = None,
    radii: float | None = None,
    rotation: tuple[int, int, int] = (0, 0, 0),
) -> tuple[Figure, Axes]:
    fig, ax = get_figure(ax)
    plot_atoms_ase(
        atoms,
        ax,
        radii=radii,
        rotation=(f"{rotation[0]}x,{rotation[1]}y,{rotation[2]}z"),
    )

    # Set axis labels with units
    ax.set_xlabel("Å")  # type: ignore bad library
    ax.set_ylabel("Å")  # type: ignore bad library

    return fig, ax


def show_legend(ax: Axes) -> None:
    ax.legend()  # type: ignore bad library


def save_fig(fig: Figure, save_dir: Path) -> None:
    fig.savefig(save_dir)  # type: ignore bad library
