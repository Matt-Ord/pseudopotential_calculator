from pathlib import Path, PosixPath

from pseudopotential_calculator.calculations.adsorbate import (
    plot_energy_against_position,
)
from pseudopotential_calculator.castep import load_all_calculators
from pseudopotential_calculator.util import save_fig

SAVE_DIR = Path(__file__).parent / "figures"

ENERGY_MAP_PATH = PosixPath("data/copper/bulk/cutoff_energy")


def _analyze_surface_energy(data_path: Path) -> None:
    calculators = load_all_calculators(data_path)
    fig = plot_energy_against_position(calculators)
    plot_name = "surface_energy.png"
    save_fig(fig, SAVE_DIR / plot_name)


if __name__ == "__main__":
    _analyze_surface_energy(ENERGY_MAP_PATH)
