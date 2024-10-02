from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atom, Atoms


def repeat_cell(cell: Atoms, amount: tuple[int, int, int]) -> Atoms:
    return cell.repeat(amount)  # type: ignore bad lib


Vector = tuple[float, float, float]


def get_basis_vectors(atoms: Atoms) -> tuple[Vector, Vector, Vector]:
    x = Vector(atoms.cell.array[0])
    y = Vector(atoms.cell.array[1])
    z = Vector(atoms.cell.array[2])
    return (x, y, z)


def add_vectors(
    vector_1: tuple[float, float, float],
    vector_2: tuple[float, float, float],
) -> tuple[float, float, float]:
    return (
        vector_1[0] + vector_2[0],
        vector_1[1] + vector_2[1],
        vector_1[2] + vector_2[2],
    )


def append_atom_to_atoms(atom: Atom, atoms: Atoms) -> Atoms:
    atoms_with_atom = atoms.copy()
    atoms_with_atom.append(atom)  # type: ignore bad lib)
    return atoms_with_atom


@dataclass
class Adsorbate:
    """class for adsorbate parameters."""

    name: str
    symbol: str
    height: float
    cut_off_energy: int
