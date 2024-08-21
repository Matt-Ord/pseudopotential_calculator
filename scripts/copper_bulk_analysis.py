import matplotlib.pyplot as plt
import numpy as np
from ase.io import read


# Test Convergence
def convergence_test() -> None:
    for n in range(1, 11):
        atoms = read(
            f"./data/bulk_cu_{n}x{n}x{n}/bulk_cu_{n}x{n}x{n}-000000.castep",
        )
        energy = atoms.get_potential_energy()

        previous_energy = None
        convergence_threshold = 1e-5  # e.g., 10^-5 eV
        is_converged = False

        if previous_energy is not None:
            energy_diff = abs(energy - previous_energy)
            print(f"k = {k}, Energy = {energy:.6f} eV, ΔE = {energy_diff:.6e} eV")

            if energy_diff < convergence_threshold:
                is_converged = True
                print(f"Converged at k = {k}")
                break
            else:
                print(f"k = {k}, Energy = {energy:.6f} eV")
        previous_energy = energy

        if calculation._error:
            print("Found error in input")
            print(calculation._error)

        if not is_converged:
            print("The energy did not converge within the tested k-point grids.")


# Plot Potential Energy against k_points
def potential_energy_plot() -> None:
    potential_energies = []
    k_points = []
    # Load the CASTEP output file
    for n in range(1, 11):
        atoms = read(
            f"./data/bulk_cu_{n}x{n}x{n}/bulk_cu_{n}x{n}x{n}-000000.castep",
        )
        potential_energies.append(atoms.get_potential_energy())
        k_points.append(n)
    print(potential_energies)
    print(k_points)
    # plot
    plt.plot(k_points, potential_energies, "o")

    # Add labels and title
    plt.xlabel("k points")
    plt.ylabel("Potential Energy")
    plt.title("Pptential Energy v.s. k points")

    # Show the plot
    plt.show()
    (plt.savefig("Potential Energy v.s. k points.png"),)

    # Plot bond length against k_points


def bond_length_plot() -> None:
    # Calculate bond length

    def vector_magnitude(vector):
        # Convert list to numpy array
        np_vector = np.array(vector)
        # Calculate the magnitude
        return np.linalg.norm(np_vector)

    bond_length_list = []
    for n in range(1, 11):
        atoms = read(
            f"/workspaces/pseudopotential_calculator/data/bulk_cu_{n}x{n}x{n}/bulk_cu_{n}x{n}x{n}-000000.cell",
        )

        bond_length_list.append(vector_magnitude(atoms.get_cell()))

    print(bond_length_list)
    # plot
    plt.plot(k_points, bond_length_list, "o")

    # Add labels and title
    plt.xlabel("Bond length")
    plt.ylabel("k points")
    plt.title("Basic Line Plot")
    plt.ylim(4.42744, 4.42746)

    # Show the plot
    plt.show()
    plt.savefig("Bond length v.s. k points.png")

    print(atoms.get_potential_energy())
