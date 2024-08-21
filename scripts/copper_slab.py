from ase.build import fcc111
from ase.constraints import FixAtoms

# Parameters
m = 5  # Number of atomic layers in the z-direction (height)
n = 3  # Number of atoms along the x and y directions (width)
L = 10  # Number of vacuum layers in the z-direction

# Create a slab of copper
# fcc111 builds an FCC (111) surface, Cu is FCC and has a lattice constant of 3.615 Å
slab = fcc111("Cu", size=(n, n, m), vacuum=L, a=3.615)

# Optionally fix the bottom layers of the slab (if needed for simulation stability)
# Uncomment the following lines to fix the bottom layer
constraint = FixAtoms(mask=[atom.tag < m - 1 for atom in slab])
slab.set_constraint(constraint)

# Optionally save the slab to a file
slab.write("Cu_slab.xyz")

fdsafds
