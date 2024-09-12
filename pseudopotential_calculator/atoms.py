from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


def repeat_cell(cell: Atoms, amount: tuple[int, int, int]) -> Atoms:
    return cell.repeat(amount)  # type: ignore bad lib


Vector = tuple[float, float, float]


def get_basis_vectors(atoms: Atoms) -> tuple[Vector, Vector, Vector]:
    x = Vector(atoms.cell.array[0])
    y = Vector(atoms.cell.array[1])
    z = Vector(atoms.cell.array[2])
    return (x, y, z)
