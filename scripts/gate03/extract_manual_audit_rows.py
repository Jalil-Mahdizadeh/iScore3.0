#!/usr/bin/env python3
"""Trace deterministic Gate-3 audit measurements to raw BindingDB rows."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import immutable_write, sha256_file, stable_json_bytes  # noqa: E402
from iscore3.gate03.curation import serialize_tsv  # noqa: E402


RAW_FIELDS = (
    "BindingDB Reactant_set_id",
    "BindingDB Ligand Name",
    "Target Name",
    "Kd (nM)",
    "pH",
    "Temp (C)",
    "Curation/DataSource",
    "Article DOI",
    "PMID",
    "Authors",
    "Institution",
    "Date of publication",
    "PDB ID(s) for Ligand-Target Complex",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bindingdb", default="data/raw/bindingdb/202608/BindingDB_All.tsv"
    )
    parser.add_argument(
        "--dataset", default="data/processed/gate03/gate03-strict-kd-v2.tsv"
    )
    parser.add_argument(
        "--config", default="configs/gate03/manual-primary-source-audit-v4.yaml"
    )
    parser.add_argument("--version", default="v4")
    args = parser.parse_args()
    bindingdb = (ROOT / args.bindingdb).resolve()
    dataset_path = (ROOT / args.dataset).resolve()
    config_path = (ROOT / args.config).resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    selections = {row["series_id"]: row for row in config["selected_series"]}
    by_series: dict[str, list[dict[str, str]]] = {key: [] for key in selections}
    for row in read_tsv(dataset_path):
        if row["series_id"] in by_series:
            by_series[row["series_id"]].append(row)
    missing = sorted(key for key, rows in by_series.items() if not rows)
    if missing:
        raise RuntimeError(f"Configured audit series missing from dataset: {missing}")

    selected_observations: list[tuple[str, dict[str, str]]] = []
    source_ids: set[str] = set()
    for series_id in sorted(by_series):
        rows = by_series[series_id]
        low = min(rows, key=lambda row: (float(row["pKd"]), row["observation_id"]))
        high = max(rows, key=lambda row: (float(row["pKd"]), row["observation_id"]))
        selected_observations.extend((("lowest_pKd", low), ("highest_pKd", high)))
        for row in (low, high):
            source_ids.update(row["source_reactant_set_ids"].split(";"))

    raw_rows: dict[str, dict[str, str]] = {}
    with bindingdb.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_fields = set(RAW_FIELDS).difference(reader.fieldnames or ())
        if missing_fields:
            raise RuntimeError(f"BindingDB fields missing: {sorted(missing_fields)}")
        for row in reader:
            reactant_id = row["BindingDB Reactant_set_id"].strip()
            if reactant_id in source_ids:
                raw_rows[reactant_id] = row
                if len(raw_rows) == len(source_ids):
                    break
    unmatched = sorted(source_ids.difference(raw_rows))
    if unmatched:
        raise RuntimeError(f"Raw BindingDB source rows not found: {unmatched}")

    output: list[dict[str, str]] = []
    for extreme, observation in selected_observations:
        selection = selections[observation["series_id"]]
        ids = observation["source_reactant_set_ids"].split(";")
        for source_index, reactant_id in enumerate(ids, start=1):
            raw = raw_rows[reactant_id]
            output.append(
                {
                    "audit_id": f"G3AUD-{len(output) + 1:03d}",
                    "series_id": observation["series_id"],
                    "observation_id": observation["observation_id"],
                    "affinity_extreme": extreme,
                    "source_replicate_index": str(source_index),
                    "consensus_pKd": observation["pKd"],
                    "consensus_source_Kd_nM": observation["source_Kd_nM"],
                    "ligand_inchikey": observation["ligand_id"],
                    "canonical_smiles": observation["canonical_smiles"],
                    "publication_id": observation["publication_id"],
                    "target_class": selection["target_class"],
                    "publication_era": selection["publication_era"],
                    "site_reference_pdb_id": observation["site_reference_pdb_id"],
                    "site_reference_ligand_comp_id": observation[
                        "site_reference_ligand_comp_id"
                    ],
                    **{f"bindingdb_{key}": raw[key] for key in RAW_FIELDS},
                }
            )
    output.sort(key=lambda row: row["audit_id"])
    output_path = ROOT / (
        f"reports/gate03/evidence/primary-source-audit-candidates-{args.version}.tsv"
    )
    immutable_write(output_path, serialize_tsv(output))
    audit = {
        "schema_version": 1,
        "status": "candidate_rows_extracted_not_yet_primary_source_verified",
        "selection_config": {"path": str(config_path), "sha256": sha256_file(config_path)},
        "dataset": {"path": str(dataset_path), "sha256": sha256_file(dataset_path)},
        "bindingdb": {"path": str(bindingdb), "sha256": sha256_file(bindingdb)},
        "selected_publications": len(selections),
        "selected_consensus_measurements": len(selected_observations),
        "raw_source_rows": len(output),
        "raw_source_ids_complete": not unmatched,
    }
    immutable_write(
        ROOT
        / f"reports/gate03/evidence/primary-source-audit-candidates-{args.version}.json",
        stable_json_bytes(audit),
    )
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
