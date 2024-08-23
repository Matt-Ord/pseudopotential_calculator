import config
from ase.io import read


# Test Convergence
def convergence_test() -> None:
    for n in config.k_range:
        atoms = read(
            f"./data/bulk_cu_{n}x{n}x{n}/bulk_cu_{n}x{n}x{n}.castep",
        )
        energy = atoms.get_cell_lengths_and_angles()
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

.get_cell_lengths_and_angles()
