from ase.build import surface
from ase.io import read, write

if __name__ == "__main__":
    bulk_copper = read("/Users/tychus/summer_project/Cu.cell", format="castep-cell")
    slab_copper = surface(bulk_copper, (0, 0, 1), 6, vacuum=3)
    # write the model file to be viewed
    write("/Users/tychus/summer_project/Cu_out.cell", b, format="castep-cell")

    _prepare_vaccum_layer_convergence(slab_copper, K_POINTS_PATH_WC)
