from pathlib import Path

from pseudopotential_calculator.calculations.bulk import (
    plot_cell_length_convergence,
    plot_energy_convergence,
)
from pseudopotential_calculator.castep import load_all_calculators

if __name__ == "__main__":
    data_path = Path("data/copper/bulk")

    calculators = load_all_calculators(data_path)
    fig, _, _ = plot_cell_length_convergence(calculators)
    fig.show()

    fig, _, _ = plot_energy_convergence(calculators)
    fig.tight_layout()
    fig.show()
    input()
