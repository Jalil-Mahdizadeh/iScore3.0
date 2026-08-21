#!/usr/bin/env python3
"""Finalize all 72 Davis chemical-identity dispositions, without labels."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import xlrd

from iscore3.gate4a.admission import (
    ALLOWED_PARENT_SALT_ROWS,
    AMBIGUOUS_DERIVATIVE_ROWS,
    DAVIS_PDF_PAGE_BY_SOURCE_ROW,
    EXPLICIT_SOURCE_STEREO_ROWS,
    UNRESOLVED_STEREO_ROWS,
    inchikey_connectivity,
)
from iscore3.provenance import load_json, verify_source_manifest


def _matches(record: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        match
        for query in record.get("queries", [])
        for match in query.get("matches", [])
        if isinstance(match, dict)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--candidate-manifest",
        type=Path,
        default=Path("data/manifests/gate4a/pubchem-candidates-v1.json"),
    )
    parser.add_argument(
        "--candidate-evidence",
        type=Path,
        default=Path("data/raw/gate4a/pubchem/davis-compound-candidates-v1.json"),
    )
    parser.add_argument(
        "--affinity-table",
        type=Path,
        default=Path("data/raw/gate4a/davis2011/supplementary_table_4.xls"),
    )
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=Path("data/manifests/gate4a/source-files-v1.json"),
    )
    parser.add_argument(
        "--chembl-manifest",
        type=Path,
        default=Path("data/manifests/gate4a/chembl-identity-evidence-v1.json"),
    )
    parser.add_argument(
        "--chembl-evidence",
        type=Path,
        default=Path("data/raw/gate4a/chembl/davis-identity-evidence-v1.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/gate4a/davis-compound-adjudication-v1.tsv"),
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=Path("reports/gate4a/evidence/davis-compound-admission-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verify_source_manifest(args.candidate_manifest, repository_root=root)
    verify_source_manifest(args.chembl_manifest, repository_root=root)
    verify_source_manifest(args.source_manifest, repository_root=root)
    candidate_payload = load_json(root / args.candidate_evidence)
    candidate_records = sorted(candidate_payload["records"], key=lambda row: row["source_row"])
    matrix = xlrd.open_workbook(
        str(root / args.affinity_table), on_demand=True
    ).sheet_by_index(0)
    rows: list[dict[str, Any]] = []
    for record in candidate_records:
        source_row = int(record["source_row"])
        properties = record.get("candidate_properties", [])
        candidate = properties[0] if len(properties) == 1 else {}
        rows.append(
            {
                "source_row": source_row,
                "source_name": record["source_name"],
                "affinity_matrix_name": str(matrix.cell_value(0, source_row + 1)),
                "pubchem_cid": candidate.get("CID", ""),
                "pubchem_inchikey": candidate.get("InChIKey", ""),
                "pubchem_smiles": candidate.get("SMILES", ""),
            }
        )
    evidence = load_json(root / args.chembl_evidence)
    by_row = {int(record["source_row"]): record for record in evidence["records"]}
    if len(rows) != 72 or len(by_row) != 72:
        raise RuntimeError("compound adjudication requires 72 ordered Davis records")

    output_rows: list[dict[str, Any]] = []
    for row in rows:
        source_row = int(row["source_row"])
        base = {
            "source_row": source_row,
            "source_name": row["source_name"],
            "affinity_matrix_name": row["affinity_matrix_name"],
            "publication_doi": "10.1038/nbt.1990",
            "publication_structure_pdf_page": DAVIS_PDF_PAGE_BY_SOURCE_ROW[source_row],
            "pubchem_cid": row["pubchem_cid"],
            "pubchem_inchikey": row["pubchem_inchikey"],
            "reviewer": "OpenAI Codex (committee verification required before label release)",
            "review_method": "manual side-by-side 2D connectivity/stereo review plus registry keys",
        }
        if source_row in AMBIGUOUS_DERIVATIVE_ROWS:
            output_rows.append(
                {
                    **base,
                    "decision": "QUARANTINED",
                    "decision_reason": "unspecified derivative; no unique molecular identity",
                    "publication_connectivity_status": "unresolvable",
                    "stereochemistry_status": "unresolvable",
                    "salt_status": "unresolvable",
                    "protomer_policy": "not_applicable",
                    "model_parent_smiles": "",
                    "model_parent_inchikey": "",
                    "chembl_agreement": "not_available",
                    "chembl_ids": "",
                    "independent_identity_evidence": "",
                }
            )
            continue

        molecule = Chem.MolFromSmiles(row["pubchem_smiles"])
        if molecule is None:
            raise RuntimeError(f"invalid PubChem candidate at source row {source_row}")
        fragment_count = len(Chem.GetMolFrags(molecule))
        if fragment_count > 1 and source_row not in ALLOWED_PARENT_SALT_ROWS:
            raise RuntimeError(f"unexpected multicomponent record at source row {source_row}")
        parent = rdMolStandardize.FragmentParent(molecule)
        parent_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        parent_inchikey = Chem.MolToInchiKey(parent)
        potential_stereo = list(Chem.FindPotentialStereo(parent))
        unspecified = sum(str(item.specified) == "Unspecified" for item in potential_stereo)
        if source_row in UNRESOLVED_STEREO_ROWS and unspecified == 0:
            raise RuntimeError(f"expected unresolved stereo at source row {source_row}")

        matches = _matches(by_row[source_row])
        exact = [m for m in matches if m.get("standard_inchi_key") == parent_inchikey]
        connectivity = [
            m
            for m in matches
            if m.get("standard_inchi_key")
            and inchikey_connectivity(str(m["standard_inchi_key"]))
            == inchikey_connectivity(parent_inchikey)
        ]
        if exact:
            chembl_status = "exact_stereo_parent"
        elif connectivity:
            chembl_status = "connectivity_parent"
        else:
            chembl_status = "no_exact_synonym_record"
        chembl_records = exact or connectivity
        chembl_ids = sorted({str(record["molecule_chembl_id"]) for record in chembl_records})

        if source_row in UNRESOLVED_STEREO_ROWS:
            decision = "QUARANTINED"
            reason = (
                "publication and registries leave a potential stereogenic element unspecified"
            )
            stereo_status = "unresolved_potential_stereo"
        else:
            decision = "ACCEPTED_PARENT"
            reason = (
                "publication connectivity matches PubChem parent; "
                "stereochemistry disposition resolved"
            )
            if source_row == 2:
                stereo_status = "source_not_drawn_external_absolute_registry_defined"
            elif source_row in EXPLICIT_SOURCE_STEREO_ROWS:
                stereo_status = "publication_drawn_matches_registry"
            else:
                stereo_status = "no_unresolved_stereogenic_element"
        if source_row == 40:
            reason += (
                "; INCB018424/INCB18424 typo resolved by row position and "
                "ruxolitinib identity"
            )

        evidence_urls = [
            f"https://pubchem.ncbi.nlm.nih.gov/compound/{row['pubchem_cid']}"
        ]
        evidence_urls.extend(
            f"https://www.ebi.ac.uk/chembl/explore/compound/{identifier}"
            for identifier in chembl_ids
        )
        if source_row == 2:
            evidence_urls.append("https://drugs.ncats.io/drug/3W2X0WGW6C")
        output_rows.append(
            {
                **base,
                "decision": decision,
                "decision_reason": reason,
                "publication_connectivity_status": (
                    "consistent_parent_counterion_not_drawn" if fragment_count > 1 else "consistent"
                ),
                "stereochemistry_status": stereo_status,
                "salt_status": (
                    "hydrochloride_counterion_removed_to_publication_parent"
                    if fragment_count > 1
                    else "single_covalent_component"
                ),
                "protomer_policy": (
                    "registry_parent_as_drawn; no pH/protomer enumeration in primary"
                ),
                "model_parent_smiles": parent_smiles if decision.startswith("ACCEPTED") else "",
                "model_parent_inchikey": parent_inchikey if decision.startswith("ACCEPTED") else "",
                "chembl_agreement": chembl_status,
                "chembl_ids": ";".join(chembl_ids),
                "independent_identity_evidence": ";".join(evidence_urls),
            }
        )

    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite adjudication: {output}")
    columns = list(output_rows[0])
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    decision_counts = Counter(row["decision"] for row in output_rows)
    audit = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "manual_rows_adjudicated": len(output_rows),
        "decision_counts": dict(sorted(decision_counts.items())),
        "publication_connectivity_reviewed_count": sum(
            row["publication_connectivity_status"].startswith("consistent") for row in output_rows
        ),
        "explicit_publication_stereo_match_count": sum(
            row["stereochemistry_status"] == "publication_drawn_matches_registry"
            for row in output_rows
        ),
        "unresolved_stereo_quarantine_count": sum(
            row["stereochemistry_status"] == "unresolved_potential_stereo" for row in output_rows
        ),
        "salt_parent_transform_count": sum(
            row["salt_status"] == "hydrochloride_counterion_removed_to_publication_parent"
            for row in output_rows
        ),
        "chembl_agreement_counts": dict(
            sorted(Counter(row["chembl_agreement"] for row in output_rows).items())
        ),
        "accepted_parent_count": decision_counts["ACCEPTED_PARENT"],
        "quarantined_count": decision_counts["QUARANTINED"],
        "all_72_received_explicit_disposition": len(output_rows) == 72,
        "human_committee_secondary_review_required": True,
        "information_boundary": {
            "affinity_labels_accessed": False,
            "protein_or_pocket_data_accessed": False,
            "crystallographic_ligand_coordinates_accessed": False,
        },
    }
    audit_output = root / args.audit_output
    audit_output.parent.mkdir(parents=True, exist_ok=True)
    if audit_output.exists():
        raise FileExistsError(f"refusing to overwrite adjudication audit: {audit_output}")
    audit_output.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"Davis compounds: {decision_counts['ACCEPTED_PARENT']} accepted / "
        f"{decision_counts['QUARANTINED']} quarantined"
    )


if __name__ == "__main__":
    main()
