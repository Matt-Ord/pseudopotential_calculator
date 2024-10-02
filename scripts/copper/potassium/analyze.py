from pathlib import Path, PosixPath

from pseudopotential_calculator.atoms import Adsorbate
from pseudopotential_calculator.calculations.slab import (
    plot_energy_against_slab_width,
    plot_energy_against_vacuum_height,
)
from pseudopotential_calculator.castep import (
    load_all_calculators,
)
from pseudopotential_calculator.scripting import maybe_copy_files_from_hpc
from pseudopotential_calculator.util import save_fig

adsorbate_config = Adsorbate(
    name="potassium",
    symbol="K",
    height=2.0,
    cut_off_energy=600,
)

VACUUM_HEIGHT_PATH = PosixPath(f"data/copper/{adsorbate_config.name}/vacuum_height")
SLAB_WIDTH_PATH = PosixPath(f"data/copper/{adsorbate_config.name}/width")
ENERGY_MAP_PATH = PosixPath(f"data/copper/{adsorbate_config.name}/energy_map")
SAVE_DIR = Path(__file__).parent / "figures"


def _analyze_convergence_with_vacuum_height(
    data_path: PosixPath,
) -> None:
    maybe_copy_files_from_hpc(data_path, data_path)
    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_vacuum_height(
        calculators,
    )
    fig.tight_layout()
    plot_name = "energy_against_vacuum_height.png"
    save_fig(fig, SAVE_DIR / plot_name)


def _analyze_convergence_with_slab_width(data_path: PosixPath) -> None:
    maybe_copy_files_from_hpc(data_path, data_path)
    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_slab_width(calculators)
    fig.tight_layout()
    plot_name = "energy_against_slab_width.png"
    save_fig(fig, SAVE_DIR / plot_name)


def _get_energy_map_from_hpc(data_path: PosixPath) -> None:
    maybe_copy_files_from_hpc(data_path, data_path)


if __name__ == "__main__":
    _analyze_convergence_with_slab_width(SLAB_WIDTH_PATH)
    _analyze_convergence_with_vacuum_height(VACUUM_HEIGHT_PATH)
    _get_energy_map_from_hpc(ENERGY_MAP_PATH)
