#!/usr/bin/env python3
"""Build selected Gate-3 CA views and run pinned US-align leakage checks."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import argparse

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import (  # noqa: E402
    immutable_write,
    sha256_file,
    stable_json_bytes,
    utc_now,
)
from iscore3.gate03.structure_mapping import read_tsv  # noqa: E402
from iscore3.protein.pocket_features import read_mmcif_atoms  # noqa: E402
from iscore3.protein.structure_views import (  # noqa: E402
    build_structural_similarity,
    write_ca_pdb,
)


def build_views(version: str) -> Path:
    selection_path = (
        ROOT
        / f"data/splits/gate03/prestructure-independent-selection-{version}.tsv"
    )
    mapping_path = (
        ROOT / "data/processed/gate03/strict-structure-site-mappings-v2.tsv"
    )
    sites_path = ROOT / "data/manifests/gate03-reference-sites-v2.json"
    selected = {row["series_id"] for row in read_tsv(selection_path)}
    mappings = {
        row["series_id"]: row
        for row in read_tsv(mapping_path)
        if row["series_id"] in selected
    }
    sites = {
        row["construct_group_id"]: row
        for row in json.loads(sites_path.read_text(encoding="utf-8"))["definitions"]
        if row["construct_group_id"] in selected
    }
    if set(mappings) != selected or set(sites) != selected:
        raise ValueError("Selected series are missing strict mapping/site definitions")
    derived_root = ROOT / f"data/interim/gate03/usalign-views-{version}"
    views = []
    for series_id in sorted(selected):
        mapping = mappings[series_id]
        site = sites[series_id]
        structure_path = Path(mapping["structure_path"])
        atoms = read_mmcif_atoms(structure_path)
        entity_id = mapping["protein_entity_id"]
        asym_id = site["reference_protein_asym_id"]
        chain_atoms = [
            atom
            for atom in atoms
            if atom.entity_id == entity_id and atom.asym_id == asym_id
        ]
        positions = {int(value) for value in site["positions_label_seq_id"]}
        pocket_atoms = [
            atom for atom in chain_atoms if atom.seq_id in positions
        ]
        global_view = write_ca_pdb(
            chain_atoms, derived_root / f"{series_id}-S1-global.pdb"
        )
        pocket_view = write_ca_pdb(
            pocket_atoms, derived_root / f"{series_id}-S1-pocket.pdb"
        )
        views.append(
            {
                "construct_group_id": series_id,
                "view": "S1",
                "status": "eligible",
                "source_id": mapping["pdb_id"],
                "source_structure_path": str(structure_path),
                "source_structure_sha256": mapping["structure_sha256"],
                "protein_entity_id": entity_id,
                "protein_asym_id": asym_id,
                "global": global_view,
                "pocket": pocket_view,
                "site_positions": sorted(positions),
            }
        )
    manifest_path = ROOT / f"data/manifests/gate03-usalign-views-{version}.json"
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "information_boundary": {
            "query_ligand_coordinates_read": False,
            "site_source": "frozen historical-reference site only",
            "purpose": "structural leakage graph only",
        },
        "inputs": {
            "selection": {
                "path": str(selection_path),
                "sha256": sha256_file(selection_path),
            },
            "mappings": {
                "path": str(mapping_path),
                "sha256": sha256_file(mapping_path),
            },
            "sites": {"path": str(sites_path), "sha256": sha256_file(sites_path)},
        },
        "counts": {"S1_eligible_groups": len(views)},
        "views": views,
    }
    immutable_write(manifest_path, stable_json_bytes(manifest))
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default="v2")
    args = parser.parse_args()
    view_manifest = build_views(args.version)
    gate2_config = yaml.safe_load(
        (ROOT / "configs/gate02/feasibility-effective-v3.yaml").read_text(
            encoding="utf-8"
        )
    )
    result = build_structural_similarity(
        view_manifest,
        gate2_config,
        ROOT / "third_party/source_cache/usalign/USalign",
        ROOT
        / f"data/processed/gate03/structural-similarity-allpairs-{args.version}.tsv",
        ROOT
        / f"data/processed/gate03/structural-leakage-edges-{args.version}.tsv",
        ROOT
        / f"reports/gate03/evidence/structural-similarity-audit-{args.version}.json",
        workers=32,
    )
    print(json.dumps(result["counts"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
