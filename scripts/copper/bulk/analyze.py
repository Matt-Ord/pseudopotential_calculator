from pathlib import Path

from pseudopotential_calculator.calculations.bulk import plot_bond_length_convergence
from pseudopotential_calculator.castep import load_all_calculators

if __name__ == "__main__":
    data_path = Path("data")

    calculators = load_all_calculators(data_path)
    fig, _, _ = plot_bond_length_convergence(calculators)
    fig.show()
