#!/usr/bin/env python3
"""Close the admitted Davis identities under explicit project-owner authorization."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from rdkit import Chem


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path("data/processed/gate4a/davis-compound-adjudication-v1.tsv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/gate4a/davis-compound-identity-final-v2.tsv"))
    parser.add_argument("--manifest", type=Path, default=Path("reports/gate4a/evidence/davis-ligand-identity-freeze-v2.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    source_path, output_path = root / args.source, root / args.output
    with source_path.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    admitted = [row for row in source_rows if row["decision"] == "ACCEPTED_PARENT"]
    if len(admitted) != 69:
        raise RuntimeError(f"expected 69 admitted ligands, observed {len(admitted)}")

    output_rows = []
    for row in admitted:
        molecule = Chem.MolFromSmiles(row["model_parent_smiles"])
        if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
            raise RuntimeError(f"invalid final parent for source row {row['source_row']}")
        canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        inchikey = Chem.MolToInchiKey(molecule)
        if canonical != row["model_parent_smiles"] or inchikey != row["model_parent_inchikey"]:
            raise RuntimeError(f"identity round-trip failed for source row {row['source_row']}")
        output_rows.append({
            **row,
            "qa_closure_status": "FINAL_OWNER_AUTHORIZED",
            "qa_closure_basis": "explicit_project_owner_instruction_2026-08-21",
            "additional_external_manual_signoff_required": "NO",
            "historical_secondary_qa_packet_disposition": "SUPERSEDED_RETAINED_FOR_PROVENANCE",
        })

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    manifest = {
        "schema_version": 2,
        "phase": "gate4a_delta3d_ligand",
        "decision": "PASS_FINAL",
        "admitted_ligands": len(output_rows),
        "technical_validation": "PASS_69_OF_69",
        "closure_basis": "The project owner explicitly accepted the completed adjudication and technical validation as final and removed the prior external/manual sign-off requirement.",
        "additional_external_manual_chemical_signoff_required": False,
        "historical_v1_ledger": {"path": str(args.source), "sha256": sha256(source_path), "disposition": "immutable_provenance"},
        "historical_secondary_packet": {"path": "data/processed/gate4a/davis-ligand-secondary-qa-packet-v1.tsv", "disposition": "superseded_not_deleted"},
        "final_v2_ledger": {"path": str(args.output), "sha256": sha256(output_path), "rows": len(output_rows)},
        "information_boundary": {"affinity_labels_accessed": False, "receptor_data_accessed": False, "bound_ligand_coordinates_accessed": False},
    }
    manifest_path = root / args.manifest
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite {manifest_path}")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
