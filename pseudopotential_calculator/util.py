from __future__ import annotations

from typing import TYPE_CHECKING, cast

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

if TYPE_CHECKING:
    from pathlib import Path


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
