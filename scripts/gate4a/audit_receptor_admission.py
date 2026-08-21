#!/usr/bin/env python3
"""Apply the frozen standardized WT reference-domain receptor estimand."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from typing import Any

import xlrd

from iscore3.gate4a.admission import (
    ordered_identity_fraction,
    receptor_exclusion_reasons,
    select_canonical_kinase_domain,
)
from iscore3.provenance import load_json, verify_source_manifest


def _domain_sequence(record: dict[str, Any], feature: dict[str, Any] | None) -> str:
    if feature is None:
        return ""
    sequence_payload = record.get("sequence") or {}
    sequence = str(sequence_payload.get("sequence", ""))
    try:
        begin = int(feature["begin"])
        end = int(feature["end"])
    except (KeyError, TypeError, ValueError):
        return ""
    return sequence[begin - 1 : end]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--target-table",
        type=Path,
        default=Path("data/raw/gate4a/davis2011/supplementary_table_1.xls"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/receptor-reference-evidence-v1.json"),
    )
    parser.add_argument(
        "--evidence",
        type=Path,
        default=Path("data/raw/gate4a/receptors/reference-evidence-v1.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/gate4a/davis-receptor-admission-v1.tsv"),
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=Path("reports/gate4a/evidence/davis-receptor-admission-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verify_source_manifest(args.manifest, repository_root=root)
    evidence = load_json(root / args.evidence)

    names_by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in evidence["klifs_names"]:
        names_by_gene[str(record["gene_name"])].append(record)
    info_by_id = {int(record["kinase_ID"]): record for record in evidence["klifs_information"]}
    uniprot_by_accession = {
        str(record["requested_accession"]): record for record in evidence["uniprot_records"]
    }
    alphafold_by_accession = {
        str(record["accession"]): record for record in evidence["alphafold_records"]
    }

    sheet = xlrd.open_workbook(str(root / args.target_table), on_demand=True).sheet_by_index(0)
    output_rows: list[dict[str, Any]] = []
    for matrix_row in range(1, sheet.nrows):
        source_accession, gene, label, mutant, group = [
            str(sheet.cell_value(matrix_row, column)).strip() for column in range(5)
        ]
        matches = names_by_gene.get(gene, [])
        match = matches[0] if len(matches) == 1 else {}
        kinase_id = int(match["kinase_ID"]) if match else None
        info = info_by_id.get(kinase_id, {}) if kinase_id is not None else {}
        uniprot_accession = str(match.get("accession", ""))
        uniprot = uniprot_by_accession.get(uniprot_accession, {})
        alphafold = alphafold_by_accession.get(uniprot_accession, {})
        domain_status, domain = select_canonical_kinase_domain(
            label, uniprot.get("domain_features", [])
        )
        pocket = str(info.get("pocket", "")) if isinstance(info.get("pocket", ""), str) else ""
        selected_domain_sequence = _domain_sequence(uniprot, dict(domain) if domain else None)
        selected_score = ordered_identity_fraction(pocket, selected_domain_sequence)
        all_core_domains = [
            feature
            for feature in uniprot.get("domain_features", [])
            if feature.get("description")
            in {"Protein kinase", "Protein kinase 1", "Protein kinase 2"}
        ]
        all_scores = [
            ordered_identity_fraction(pocket, _domain_sequence(uniprot, feature))
            for feature in all_core_domains
        ]
        best_score = max(all_scores, default=0.0)
        mapping_is_best = selected_score >= best_score - 1e-12 if domain else False

        reasons = list(
            receptor_exclusion_reasons(
                assay_label=label,
                mutant_flag=mutant,
                kinase_group=group,
                klifs_match_count=len(matches),
                pocket_length=len(pocket),
                domain_status=domain_status,
                alphafold_available=bool(alphafold.get("available", False)),
            )
        )
        if domain and (selected_score < 0.95 or not mapping_is_best):
            reasons.append("klifs_pocket_to_selected_domain_alignment_low_confidence")
        decision = "ACCEPTED_REFERENCE_DOMAIN" if not reasons else "EXCLUDED_PRIMARY"
        begin = str(domain.get("begin", "")) if domain else ""
        end = str(domain.get("end", "")) if domain else ""
        estimand_id = (
            f"UniProt:{uniprot_accession}:{begin}-{end}|KLIFS:{kinase_id}:positions1-85"
            if decision == "ACCEPTED_REFERENCE_DOMAIN"
            else ""
        )
        output_rows.append(
            {
                "matrix_row_1_based": matrix_row,
                "source_refseq_accession": source_accession,
                "gene_symbol": gene,
                "assay_target_label": label,
                "source_mutant_flag": mutant,
                "source_kinase_group": group,
                "primary_decision": decision,
                "exclusion_reasons": ";".join(dict.fromkeys(reasons)),
                "estimand_id": estimand_id,
                "klifs_kinase_id": kinase_id or "",
                "klifs_name": info.get("name", ""),
                "klifs_family": info.get("family", ""),
                "klifs_group": info.get("group", ""),
                "klifs_pocket_sequence": pocket,
                "klifs_pocket_length": len(pocket),
                "uniprot_accession": uniprot_accession,
                "uniprot_entry_id": uniprot.get("uniProtkbId", ""),
                "domain_resolution_status": domain_status,
                "canonical_domain_description": domain.get("description", "") if domain else "",
                "canonical_domain_begin": begin,
                "canonical_domain_end": end,
                "pocket_domain_ordered_identity": f"{selected_score:.6f}",
                "selected_domain_is_best_pocket_match": str(mapping_is_best).lower(),
                "alphafold_available": str(bool(alphafold.get("available", False))).lower(),
                "alphafold_entry_id": alphafold.get("entry_id", ""),
                "alphafold_version": alphafold.get("latest_version", ""),
                "alphafold_model_created_date": alphafold.get("model_created_date", ""),
                "alphafold_cif_url": alphafold.get("cif_url", ""),
                "receptor_view_role": (
                    "primary_predicted_standardized_WT_reference_domain"
                    if decision == "ACCEPTED_REFERENCE_DOMAIN"
                    else "not_admitted"
                ),
            }
        )

    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite receptor admission: {output}")
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    accepted = [row for row in output_rows if row["primary_decision"].startswith("ACCEPTED")]
    exclusion_counts = Counter(
        reason for row in output_rows for reason in row["exclusion_reasons"].split(";") if reason
    )
    audit = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "reference_estimand": (
            "reviewed human UniProt canonical WT core kinase domain; fixed KLIFS 85-position "
            "pocket; AlphaFold DB canonical-sequence prediction"
        ),
        "exact_kinomescan_construct_estimand_status": "BLOCKED_NOT_REPORTED_BY_SOURCE",
        "standardized_reference_domain_estimand_status": "PASS_FOR_ADMITTED_ROWS",
        "assay_row_count": len(output_rows),
        "accepted_assay_row_count": len(accepted),
        "excluded_assay_row_count": len(output_rows) - len(accepted),
        "accepted_unique_estimand_count": len({row["estimand_id"] for row in accepted}),
        "accepted_unique_uniprot_count": len({row["uniprot_accession"] for row in accepted}),
        "accepted_unique_klifs_pocket_count": len({row["klifs_kinase_id"] for row in accepted}),
        "exclusion_reason_counts_nonexclusive": dict(sorted(exclusion_counts.items())),
        "fixed_pocket_position_count": 85,
        "pocket_mapping_minimum_ordered_identity": 0.95,
        "accepted_pocket_alignment_minimum": min(
            (float(row["pocket_domain_ordered_identity"]) for row in accepted), default=None
        ),
        "predicted_view_coverage_accepted": sum(
            row["alphafold_available"] == "true" for row in accepted
        ),
        "apo_view_status": (
            "CANDIDATE_AUDIT_PENDING_STRICT_ATOM_LEVEL_POCKET_OCCUPANCY; predicted view is "
            "the primary non-holo view"
        ),
        "leakage_controls": {
            "query_ligand_used_to_define_pocket": False,
            "holo_coordinates_used_for_primary_view": False,
            "affinity_labels_used_for_mapping_or_exclusion": False,
            "mutants_or_state_specific_rows_collapsed_into_WT": False,
        },
    }
    audit_output = root / args.audit_output
    audit_output.parent.mkdir(parents=True, exist_ok=True)
    if audit_output.exists():
        raise FileExistsError(f"refusing to overwrite receptor audit: {audit_output}")
    audit_output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"Davis receptors: {len(accepted)}/{len(output_rows)} assay rows accepted; "
        f"{audit['accepted_unique_estimand_count']} unique reference-domain estimands"
    )


if __name__ == "__main__":
    main()
