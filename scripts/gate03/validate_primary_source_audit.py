#!/usr/bin/env python3
"""Validate the manual Gate-3 primary-source registry against frozen candidates."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import immutable_write, sha256_file, stable_json_bytes  # noqa: E402
from iscore3.gate03.curation import serialize_tsv  # noqa: E402


def main() -> None:
    registry_path = ROOT / "configs/gate03/primary-source-verification-v4.yaml"
    registry = yaml.safe_load(registry_path.read_text(encoding="utf-8"))
    candidate_path = ROOT / registry["candidate_table"]
    with candidate_path.open("r", encoding="utf-8", newline="") as handle:
        raw_candidates = list(csv.DictReader(handle, delimiter="\t"))
    candidates: dict[str, list[dict[str, str]]] = {}
    for row in raw_candidates:
        candidates.setdefault(row["observation_id"], []).append(row)
    verification = {
        observation_id: {**series, **measurement}
        for series_id, series in registry["series"].items()
        for observation_id, measurement in series["measurements"].items()
    }
    if set(candidates) != set(verification):
        raise RuntimeError(
            "Verification/candidate observation mismatch: "
            f"missing={sorted(set(candidates) - set(verification))}, "
            f"extra={sorted(set(verification) - set(candidates))}"
        )
    output = []
    for observation_id in sorted(candidates):
        candidate = candidates[observation_id][0]
        verified = verification[observation_id]
        output.append(
            {
                "observation_id": observation_id,
                "series_id": candidate["series_id"],
                "publication_id": candidate["publication_id"],
                "target_name": candidate["bindingdb_Target Name"],
                "affinity_extreme": candidate["affinity_extreme"],
                "bindingdb_Kd_nM": candidate["consensus_source_Kd_nM"],
                "ligand_identifier": verified["compound"],
                "verification_status": verified["status"],
                "primary_value": verified["primary_value"],
                "primary_location": verified["location"],
                "assay": verified["assay"],
                "target_construct": verified["target_construct"],
                "reference_structure_relation": verified[
                    "reference_structure_relation"
                ],
                "source_url": verified["source_url"],
            }
        )
    verified_rows = [row for row in output if row["verification_status"] == "verified"]
    verified_series = {row["series_id"] for row in verified_rows}
    fraction = len(verified_rows) / len(output)
    manifest_path = ROOT / registry["source_manifest"]
    checks = {
        "minimum_publications_12": len(registry["series"]) >= 12,
        "all_publications_have_at_least_one_verified_measurement": len(verified_series)
        == len(registry["series"]),
        "minimum_consensus_measurements_24": len(output) >= 24,
        "minimum_verified_fraction_0p90": fraction >= 0.90,
        "source_manifest_exists": manifest_path.is_file(),
    }
    audit = {
        "schema_version": 1,
        "status": "PASS" if all(checks.values()) else "FAIL",
        "registry": {"path": str(registry_path), "sha256": sha256_file(registry_path)},
        "candidates": {"path": str(candidate_path), "sha256": sha256_file(candidate_path)},
        "source_manifest": {
            "path": str(manifest_path),
            "sha256": sha256_file(manifest_path) if manifest_path.is_file() else "",
        },
        "publications": len(registry["series"]),
        "consensus_measurements": len(output),
        "verified_measurements": len(verified_rows),
        "verified_fraction": fraction,
        "failed_or_unresolved_measurements": len(output) - len(verified_rows),
        "checks": checks,
    }
    immutable_write(
        ROOT / "reports/gate03/evidence/primary-source-verification-v4.tsv",
        serialize_tsv(output),
    )
    immutable_write(
        ROOT / "reports/gate03/evidence/primary-source-verification-v4.json",
        stable_json_bytes(audit),
    )
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
