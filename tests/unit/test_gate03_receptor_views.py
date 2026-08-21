import numpy as np

from iscore3.gate03.receptor_views import POCKET_V2_NAMES, pocket_descriptor_v2
from iscore3.protein.pocket_features import AtomRecord, PocketInstance


def _pocket(offset=(0.0, 0.0, 0.0)):
    names = ["ALA", "ASP", "LYS", "PHE"]
    atoms = []
    for position, (name, xyz) in enumerate(
        zip(names, ((0, 0, 0), (4, 0, 0), (0, 5, 0), (0, 0, 6)), strict=True),
        start=1,
    ):
        for atom_name, shift, element in (("CA", (0, 0, 0), "C"), ("N", (0.2, 0.1, 0), "N")):
            point = tuple(float(xyz[i] + shift[i] + offset[i]) for i in range(3))
            atoms.append(
                AtomRecord(
                    "ATOM", element, atom_name, "", name, "A", "1", position,
                    point, 1.0, "A", str(position), 80.0 + position,
                )
            )
    return PocketInstance(
        "TEST", "0" * 64, "1", "A", (1, 2, 3, 4), (1, 2, 3, 4), (),
        {index + 1: name for index, name in enumerate(names)}, tuple(atoms),
    )


def test_pocket_v2_is_schema_exact_and_translation_invariant():
    first = pocket_descriptor_v2(_pocket())
    second = pocket_descriptor_v2(_pocket((100.0, -20.0, 7.5)))
    assert tuple(first) == POCKET_V2_NAMES
    assert np.allclose(list(first.values()), list(second.values()), atol=1.0e-12)
    assert np.isfinite(list(first.values())).all()
