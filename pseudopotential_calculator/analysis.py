from __future__ import annotations

import matplotlib.pyplot as plt
from ase.io import read

from pseudopotential_calculator.config import CASTEPConfig


def potential_energy_convergence_test(
    config: CASTEPConfig = CASTEPConfig(),
    convergence_threshold: float = 1e-5,
) -> tuple[list[int], list[float]]:
    k_values = []
    energies = []
    previous_energy = None
    # e.g., 10^-5 eV
    is_converged = False

    for n in config.k_range:
        atoms = read(f"./data/bulk_cu_{n}x{n}x{n}/bulk_cu_{n}x{n}x{n}.castep")
        energy = atoms.get_potential_energy()

        k_values.append(n)
        energies.append(energy)

        if previous_energy is not None:
            energy_diff = abs(energy - previous_energy)
            print(f"k = {n}, Energy = {energy:.6f} eV, ΔE = {energy_diff:.6e} eV")

            if energy_diff < convergence_threshold:
                is_converged = True
                print(f"Converged at k = {n}")
                break
        else:
            print(f"k = {n}, Energy = {energy:.6f} eV")

        previous_energy = energy

    if not is_converged:
        print("The energy did not converge within the tested k-point grids.")

    return k_values, energies


# Function to plot potential energy vs k-points
def potential_energy_plot(k_values: list[int], bond_lengths: list[float]) -> None:
    ax = plt.subplot()
    fig = ax.get_figure()
    ax.plot(k_values, bond_lengths, marker="o", linestyle="-")
    ax.set_xlabel("Potential Energy (eV)")
    ax.set_ylabel("Bond Length (Angstrom)")
    ax.set_title("Convergence Test: Potential Energy vs k-point grid")
    fig.save("potential energy vs k-points")  # Save the plot to a file
    print("Plot saved")
    plt.close()


# Function to perform convergence test and collect data on bond length
def bond_length_convergence_test(
    convergence_threshold=1e-5,
) -> tuple[list[int], list[list[float]]]:
    k_values = []
    bond_lengths = []
    previous_bond_length = None
    is_converged = False

    for n in CASTEPConfig.k_range:
        atoms = read(f"./data/bulk_cu_{n}x{n}x{n}/bulk_cu_{n}x{n}x{n}.castep")
        # Get cell parameters using deprecated method
        lengths_and_angles = atoms.get_cell_lengths_and_angles()

        k_values.append(n)
        bond_length = lengths_and_angles[0]
        print(bond_length)
        bond_lengths.append(bond_length)

        if previous_bond_length is not None:
            # Compute the difference in cell parameters
            length_diff = abs(bond_length - previous_bond_length)

            print(f"k = {n}, Bond length = {bond_length}")
            print(f"ΔLengths = {length_diff}")

            if length_diff < convergence_threshold:
                is_converged = True
                print(f"Converged at k = {n}")
                break
        else:
            print(f"k = {n}, Bond length = {bond_length}")
        previous_bond_length = bond_length

    if not is_converged:
        print("The cell parameters did not converge within the tested k-point grids.")

    return k_values, bond_lengths


# Function to plot bond length vs k-points
def bond_length_plot(k_values: list[int], bond_lengths: list[float]) -> None:
    ax = plt.subplot()
    fig = ax.get_figure()
    ax.plot(k_values, bond_lengths, marker="o", linestyle="-")
    ax.set_xlabel("k-point grid size (n)")
    ax.set_ylabel("Bond Length (Angstrom)")
    ax.set_title("Convergence Test: Bond Length vs k-point grid")
    fig.save("bond length vs k-points")  # Save the plot to a file
    print("Plot saved")
    plt.close()
