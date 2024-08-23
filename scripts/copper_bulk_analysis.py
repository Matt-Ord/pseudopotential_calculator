from pseudopotential_calculator import analysis

if __name__ == "__main__":
    k_values, energies = analysis.potential_energy_convergence_test()
    analysis.potential_energy_plot(k_values, energies)

    k_values, bond_lengths = analysis.bond_length_convergence_test()
    analysis.bond_length_plot(k_values, bond_lengths)
