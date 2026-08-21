#!/usr/bin/env python3
"""Prove that every retained BindingDB top-score tie is affinity-fact equivalent."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

from iscore3.data.bindingdb_audit import (
    _candidate_assessment,
    _pdb_ids,
    scan_bindingdb,
)
from iscore3.data.rcsb_gate01 import (
    immutable_write,
    sha256_file,
    stable_json_bytes,
    utc_now,
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict-pilot", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    strict = [row for row in rows(args.strict_pilot) if row["role"] == "supervised_s0"]
    strict_by_id = {row["observation_id"]: row for row in strict}
    provenance = {row["observation_id"]: row for row in rows(args.provenance)}
    if not set(strict_by_id).issubset(provenance):
        raise RuntimeError("Strict observations are absent from the provenance table")
    tied_ids = sorted(
        observation_id
        for observation_id in strict_by_id
        if int(provenance[observation_id]["top_score_tie_count"]) > 1
    )
    candidates, scan = scan_bindingdb(
        args.archive,
        [strict_by_id[observation_id]["pdb_id"] for observation_id in tied_ids],
    )
    by_pdb: dict[str, list[dict[str, str]]] = defaultdict(list)
    for candidate in candidates:
        for pdb_id in _pdb_ids(candidate["PDB ID(s) for Ligand-Target Complex"]):
            by_pdb[pdb_id].append(candidate)

    audits = []
    for observation_id in tied_ids:
        pilot = strict_by_id[observation_id]
        assessments = [
            _candidate_assessment(pilot, candidate)
            for candidate in by_pdb.get(pilot["pdb_id"].upper(), [])
        ]
        assessments.sort(
            key=lambda item: (
                -item["score"],
                item["candidate"]["BindingDB Reactant_set_id"],
            )
        )
        if not assessments:
            raise RuntimeError(f"No rescanned candidate for {observation_id}")
        top_score = assessments[0]["score"]
        tied = [item for item in assessments if item["score"] == top_score]
        expected_count = int(provenance[observation_id]["top_score_tie_count"])
        selected_id = provenance[observation_id]["bindingdb_reactant_set_id"]
        equivalent = bool(
            len(tied) == expected_count
            and selected_id
            in {item["candidate"]["BindingDB Reactant_set_id"] for item in tied}
            and all(
                item["high_confidence"]
                and item["pdb_match"]
                and item["ligand_match"]
                and item["uniprot_match"]
                and item["kd_match"]
                and item["measurement_publication_present"]
                for item in tied
            )
        )
        publications = sorted(
            {
                (
                    item["candidate"]["Article DOI"].strip().lower(),
                    item["candidate"]["PMID"].strip(),
                )
                for item in tied
            }
        )
        audits.append(
            {
                "observation_id": observation_id,
                "pdb_id": pilot["pdb_id"],
                "expected_tie_count": expected_count,
                "rescanned_tie_count": len(tied),
                "selected_reactant_set_id": selected_id,
                "tied_reactant_set_ids": [
                    item["candidate"]["BindingDB Reactant_set_id"] for item in tied
                ],
                "affinity_fact_equivalent": equivalent,
                "distinct_publication_count": len(publications),
                "publications": [
                    {"doi": doi, "pmid": pmid} for doi, pmid in publications
                ],
            }
        )
    passed = all(record["affinity_fact_equivalent"] for record in audits)
    report = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "PASS" if passed else "FAIL",
        "definition": (
            "all top-score alternatives independently match the same retained PDB, ligand, "
            "UniProt target, exact Kd, and a populated measurement publication"
        ),
        "inputs": {
            "strict_pilot": {
                "path": str(args.strict_pilot),
                "sha256": sha256_file(args.strict_pilot),
            },
            "provenance": {
                "path": str(args.provenance),
                "sha256": sha256_file(args.provenance),
            },
            "archive": {
                "path": str(args.archive),
                "sha256": sha256_file(args.archive),
            },
        },
        "scan": scan,
        "counts": {
            "strict_supervised_observations": len(strict),
            "retained_top_score_ties": len(tied_ids),
            "affinity_fact_equivalent_ties": sum(
                record["affinity_fact_equivalent"] for record in audits
            ),
            "ties_with_multiple_publications": sum(
                record["distinct_publication_count"] > 1 for record in audits
            ),
        },
        "records": audits,
    }
    immutable_write(args.output, stable_json_bytes(report))
    print(
        json.dumps({"status": report["status"], "counts": report["counts"]}, indent=2)
    )
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
