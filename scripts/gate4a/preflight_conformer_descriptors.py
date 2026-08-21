#!/usr/bin/env python3
"""Outcome-blind conformer descriptor preflight on unaccepted Davis candidates."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import statistics
from typing import Any

import numpy as np

from iscore3.gate4a.conformers import (
    ConformerConfig,
    ConformerGenerationError,
    generate_conformer_descriptors,
)
from iscore3.provenance import load_json, sha256_file, verify_source_manifest


def _candidate(record: dict[str, Any]) -> dict[str, Any]:
    properties = record.get("candidate_properties", [])
    return properties[0] if len(properties) == 1 and isinstance(properties[0], dict) else {}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/pubchem-candidates-v1.json"),
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=Path("data/raw/gate4a/pubchem/davis-compound-candidates-v1.json"),
    )
    parser.add_argument(
        "--feature-output",
        type=Path,
        default=Path("data/features/gate4a/davis-candidate-conformers-v1.json"),
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=Path("reports/gate4a/evidence/davis-conformer-preflight-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verified = verify_source_manifest(args.manifest, repository_root=root)
    acquisition = load_json(root / args.candidates)
    records = acquisition.get("records")
    if not isinstance(records, list) or len(records) != 72:
        raise RuntimeError("expected the frozen 72-row candidate acquisition")
    candidate_records = [
        record for record in records if record.get("mapping_state") == "candidate_requires_manual_verification"
    ]
    if len(candidate_records) != 71:
        raise RuntimeError("expected 71 unaccepted single-candidate structures")

    config = ConformerConfig()
    successes: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    schema: tuple[str, ...] | None = None
    groups: tuple[str, ...] | None = None
    for record in candidate_records:
        properties = _candidate(record)
        try:
            vector = generate_conformer_descriptors(str(properties["SMILES"]), config)
        except (ConformerGenerationError, ValueError) as exc:
            failures.append(
                {
                    "source_row": record["source_row"],
                    "source_name": record["source_name"],
                    "pubchem_cid": properties.get("CID"),
                    "error": str(exc),
                }
            )
            continue
        if schema is None:
            schema = vector.feature_names
            groups = vector.feature_groups
        if vector.feature_names != schema or vector.feature_groups != groups:
            raise RuntimeError("descriptor schema changed across molecules")
        successes.append(
            {
                "source_row": record["source_row"],
                "source_name": record["source_name"],
                "pubchem_cid": properties["CID"],
                "pubchem_inchikey": properties["InChIKey"],
                "canonical_isomeric_smiles": vector.canonical_isomeric_smiles,
                "heavy_atom_symbols": vector.heavy_atom_symbols,
                "heavy_atom_map_numbers": vector.heavy_atom_map_numbers,
                "retained_conformer_ids": vector.retained_conformer_ids,
                "retained_energies_kcal_mol": vector.retained_energies_kcal_mol,
                "boltzmann_weights": vector.boltzmann_weights,
                "generated_geometry_sha256": vector.generated_geometry_sha256,
                "unspecified_stereocentre_count": vector.unspecified_stereocentre_count,
                "values": vector.values,
            }
        )
    if not successes or schema is None or groups is None:
        raise RuntimeError("conformer preflight produced no descriptors")

    schema_payload = {"feature_names": schema, "feature_groups": groups}
    schema_sha256 = hashlib.sha256(
        json.dumps(schema_payload, separators=(",", ":"), sort_keys=True).encode("utf-8")
    ).hexdigest()
    feature_payload = {
        "schema_version": 1,
        "phase": "gate4a",
        "role": "technical preflight on unaccepted PubChem candidates; forbidden for modeling",
        "candidate_acquisition": verified[0].__dict__,
        "config": config.__dict__,
        "feature_schema_sha256": schema_sha256,
        **schema_payload,
        "successes": successes,
        "failures": failures,
    }
    feature_output = root / args.feature_output
    if feature_output.exists():
        raise FileExistsError(f"refusing to overwrite preflight features: {feature_output}")
    feature_output.parent.mkdir(parents=True, exist_ok=True)
    with feature_output.open("x", encoding="utf-8") as handle:
        json.dump(feature_payload, handle, indent=2, sort_keys=True)
        handle.write("\n")

    matrix = np.asarray([record["values"] for record in successes], dtype=np.float64)
    geometry_indices = [
        index
        for index, group in enumerate(groups)
        if group in {"shape3d", "pharmacophore3d", "conformer_diversity", "conformer_energy"}
    ]
    geometry_standard_deviation = matrix[:, geometry_indices].std(axis=0)
    retained_counts = [len(record["retained_conformer_ids"]) for record in successes]
    diversity_index = schema.index("conformer_diversity.heavy_atom_best_rmsd_mean")
    diversity = matrix[:, diversity_index]
    audit = {
        "schema_version": 1,
        "phase": "gate4a",
        "information_boundary": {
            "inputs": ["unaccepted PubChem candidate SMILES", "frozen conformer config"],
            "affinity_labels_accessed": False,
            "protein_or_pocket_data_accessed": False,
            "crystallographic_or_docked_ligand_coordinates_accessed": False,
        },
        "candidate_structures_accepted_for_modeling": False,
        "feature_file": {
            "path": str(args.feature_output),
            "bytes": feature_output.stat().st_size,
            "sha256": sha256_file(feature_output),
        },
        "feature_schema_sha256": schema_sha256,
        "feature_count": len(schema),
        "feature_group_counts": dict(sorted(Counter(groups).items())),
        "candidate_count": len(candidate_records),
        "successful_generation_count": len(successes),
        "failed_generation_count": len(failures),
        "failures": failures,
        "unique_generated_geometry_hash_count": len(
            {record["generated_geometry_sha256"] for record in successes}
        ),
        "retained_conformer_count": {
            "minimum": min(retained_counts),
            "median": statistics.median(retained_counts),
            "maximum": max(retained_counts),
            "at_least_two": sum(value >= 2 for value in retained_counts),
        },
        "molecules_with_unassigned_stereocentres": sum(
            record["unspecified_stereocentre_count"] > 0 for record in successes
        ),
        "heavy_atom_conformer_diversity": {
            "mean_of_pairwise_best_rmsd_mean": float(diversity.mean()),
            "median_of_pairwise_best_rmsd_mean": float(np.median(diversity)),
            "molecules_above_0_5_angstrom": int(np.sum(diversity > 0.5)),
        },
        "nonconstant_geometry_feature_count": int(
            np.sum(geometry_standard_deviation > 1e-8)
        ),
        "geometry_feature_count": len(geometry_indices),
        "interpretation": (
            "This checks generator coverage and representation variation only. It is not an "
            "affinity result and cannot qualify unverified compound identities."
        ),
    }
    audit_output = root / args.audit_output
    if audit_output.exists():
        raise FileExistsError(f"refusing to overwrite preflight audit: {audit_output}")
    audit_output.parent.mkdir(parents=True, exist_ok=True)
    with audit_output.open("x", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        f"Conformer preflight: {len(successes)}/{len(candidate_records)} succeeded; "
        f"{int(np.sum(geometry_standard_deviation > 1e-8))}/{len(geometry_indices)} "
        "geometry features varied"
    )


if __name__ == "__main__":
    main()
