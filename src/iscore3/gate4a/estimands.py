"""Frozen Gate-4A model/contrast registry; this module performs no fitting."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ModelSpec:
    name: str
    main_terms: tuple[str, ...]
    interaction_ranks: tuple[tuple[str, int], ...] = ()
    role: str = "primary"

    @property
    def total_interaction_rank(self) -> int:
        return sum(rank for _, rank in self.interaction_ranks)


@dataclass(frozen=True)
class EffectSpec:
    name: str
    reference_model: str
    augmented_model: str
    interpretation: str


MODEL_SPECS: dict[str, ModelSpec] = {
    "M2D": ModelSpec("M2D", ("ligand_2d",)),
    "M3D": ModelSpec("M3D", ("ligand_2d", "ligand_3d")),
    "M2D_P": ModelSpec("M2D_P", ("ligand_2d", "pocket")),
    "MA": ModelSpec("MA", ("ligand_2d", "ligand_3d", "pocket")),
    "MI2_rank8": ModelSpec(
        "MI2_rank8",
        ("ligand_2d", "ligand_3d", "pocket"),
        (("ligand_2d_x_pocket", 8),),
    ),
    "MI23_rank4_plus_rank4": ModelSpec(
        "MI23_rank4_plus_rank4",
        ("ligand_2d", "ligand_3d", "pocket"),
        (("ligand_2d_x_pocket", 4), ("ligand_3d_x_pocket", 4)),
    ),
    "MI2_rank4": ModelSpec(
        "MI2_rank4",
        ("ligand_2d", "ligand_3d", "pocket"),
        (("ligand_2d_x_pocket", 4),),
        role="sensitivity",
    ),
    "MI2_rank4_plus_I3_rank4": ModelSpec(
        "MI2_rank4_plus_I3_rank4",
        ("ligand_2d", "ligand_3d", "pocket"),
        (("ligand_2d_x_pocket", 4), ("ligand_3d_x_pocket", 4)),
        role="sensitivity",
    ),
}


EFFECT_SPECS: dict[str, EffectSpec] = {
    "delta_3d_ligand": EffectSpec(
        "delta_3d_ligand",
        "M2D",
        "M3D",
        "incremental free-ligand 3D representation value",
    ),
    "delta_pocket_additive": EffectSpec(
        "delta_pocket_additive",
        "M3D",
        "MA",
        "additive receptor calibration without ligand-pocket interaction",
    ),
    "delta_3d_x_pocket": EffectSpec(
        "delta_3d_x_pocket",
        "MI2_rank8",
        "MI23_rank4_plus_rank4",
        "capacity-matched ligand-3D-specific statistical interaction",
    ),
}


def validate_registry() -> None:
    for effect in EFFECT_SPECS.values():
        if effect.reference_model not in MODEL_SPECS or effect.augmented_model not in MODEL_SPECS:
            raise ValueError(f"effect references an unknown model: {effect.name}")
    reference = MODEL_SPECS["MI2_rank8"]
    augmented = MODEL_SPECS["MI23_rank4_plus_rank4"]
    if reference.total_interaction_rank != augmented.total_interaction_rank:
        raise ValueError("primary interaction comparison is not capacity-matched")
    if set(reference.main_terms) != set(augmented.main_terms):
        raise ValueError("primary interaction comparison has unequal main effects")
