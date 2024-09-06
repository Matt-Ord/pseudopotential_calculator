from pathlib import Path, PosixPath

from pseudopotential_calculator.calculations.slab import plot_energy_against_nxn_slab
from pseudopotential_calculator.castep import (
    load_all_calculators,
)
from pseudopotential_calculator.scripting import maybe_copy_files_from_hpc
from pseudopotential_calculator.util import save_fig

SLAB_WITH_NA_PATH = PosixPath("data/copper/slab/vacuum_layer")

SAVE_DIR = Path(__file__).parent / "figures"


def _analyze_convergence_with_nxn_slab() -> None:
    data_path = SLAB_WITH_NA_PATH
    maybe_copy_files_from_hpc(data_path, data_path)
    calculators = load_all_calculators(data_path)

    fig, _, _ = plot_energy_against_nxn_slab(calculators)
    fig.tight_layout()
    plot_name = "energy_against_nxn_slab.png"
    save_fig(fig, SAVE_DIR / plot_name)


if __name__ == "__main__":
    _analyze_convergence_with_nxn_slab()
