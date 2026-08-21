#!/usr/bin/env python3
"""Verify the admitted ligand ledger and render an external secondary-QA packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from rdkit import Chem


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", type=Path, default=Path("data/processed/gate4a/davis-compound-adjudication-v1.tsv"))
    parser.add_argument("--qa-packet", type=Path, default=Path("data/processed/gate4a/davis-ligand-secondary-qa-packet-v1.tsv"))
    parser.add_argument("--freeze", type=Path, default=Path("reports/gate4a/evidence/davis-ligand-identity-freeze-v1.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    ledger_path = root / args.ledger
    with ledger_path.open(encoding="utf-8", newline="") as handle:
        source = list(csv.DictReader(handle, delimiter="\t"))
    admitted = [row for row in source if row["decision"] == "ACCEPTED_PARENT"]
    packet: list[dict[str, Any]] = []
    failures = []
    for row in admitted:
        molecule = Chem.MolFromSmiles(row["model_parent_smiles"])
        computed_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True) if molecule else ""
        computed_key = Chem.MolToInchiKey(molecule) if molecule else ""
        exact = computed_smiles == row["model_parent_smiles"] and computed_key == row["model_parent_inchikey"]
        single_component = molecule is not None and len(Chem.GetMolFrags(molecule)) == 1
        evidence = [value for value in row["independent_identity_evidence"].split(";") if value]
        checks = {
            "rdkit_roundtrip_exact": exact,
            "single_covalent_component": single_component,
            "publication_page_present": bool(row["publication_structure_pdf_page"]),
            "pubchem_record_present": any("pubchem.ncbi.nlm.nih.gov" in value for value in evidence),
            "stereochemistry_disposition_present": bool(row["stereochemistry_status"]),
            "salt_disposition_present": bool(row["salt_status"]),
            "protomer_policy_present": bool(row["protomer_policy"]),
        }
        if not all(checks.values()):
            failures.append({"source_row": row["source_row"], "failed": [key for key, value in checks.items() if not value]})
        blind_id = "DQA-" + hashlib.sha256(f"{row['source_row']}|{row['model_parent_inchikey']}".encode()).hexdigest()[:12].upper()
        packet.append({"qa_id": blind_id, "source_row": row["source_row"], "source_name": row["source_name"], "affinity_matrix_name": row["affinity_matrix_name"], "publication_doi": row["publication_doi"], "publication_structure_pdf_page": row["publication_structure_pdf_page"], "candidate_parent_smiles": row["model_parent_smiles"], "candidate_parent_inchikey": row["model_parent_inchikey"], "pubchem_cid": row["pubchem_cid"], "chembl_ids": row["chembl_ids"], "prior_stereochemistry_disposition": row["stereochemistry_status"], "prior_salt_disposition": row["salt_status"], "prior_protomer_policy": row["protomer_policy"], "evidence_urls": row["independent_identity_evidence"], "automated_roundtrip_status": "PASS" if all(checks.values()) else "FAIL", "secondary_connectivity_decision": "PENDING", "secondary_stereochemistry_decision": "PENDING", "secondary_salt_parent_decision": "PENDING", "secondary_protomer_decision": "PENDING", "secondary_reviewer_name": "PENDING", "secondary_review_date": "PENDING", "secondary_notes": ""})
    packet_path = root / args.qa_packet
    packet_path.parent.mkdir(parents=True, exist_ok=True)
    if packet_path.exists():
        raise FileExistsError(f"refusing to overwrite {packet_path}")
    with packet_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(packet[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(packet)
    freeze = {
        "schema_version": 1,
        "phase": "gate4a_provenance_closure",
        "source_disposition_rows": len(source),
        "admitted_parent_rows": len(admitted),
        "quarantined_rows": len(source) - len(admitted),
        "technical_integrity_decision": "PASS" if not failures and len(admitted) == 69 else "FAIL",
        "technical_failures": failures,
        "independent_secondary_qa_decision": "BLOCKED_PENDING_NAMED_EXTERNAL_REVIEWER_SIGNATURE",
        "final_chemical_identity_ledger_decision": "BLOCKED",
        "rationale": "The same research agent that produced the first adjudication cannot truthfully certify its own second pass as independent. The immutable ledger hash and a row-complete review packet are frozen for a committee chemist or otherwise independent reviewer.",
        "source_ledger": {"path": str(args.ledger), "sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest()},
        "secondary_qa_packet": {"path": str(args.qa_packet), "sha256": hashlib.sha256(packet_path.read_bytes()).hexdigest(), "rows": len(packet)},
        "release_rule": "All four secondary decision columns must be PASS for all 69 rows, reviewer name/date must be non-PENDING, and the signed packet hash must be added in a new immutable manifest before ligand labels or Delta3D-ligand fitting are released.",
        "information_boundary": {"affinity_labels_accessed": False, "receptor_data_accessed": False, "crystallographic_ligand_coordinates_accessed": False},
    }
    freeze_path = root / args.freeze
    freeze_path.parent.mkdir(parents=True, exist_ok=True)
    if freeze_path.exists():
        raise FileExistsError(f"refusing to overwrite {freeze_path}")
    freeze_path.write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
