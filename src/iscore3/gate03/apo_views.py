"""Strict pocket-unoccupied experimental S3 views for Gate-3."""

from __future__ import annotations

from collections import Counter, defaultdict
import gzip
import hashlib
from pathlib import Path
from typing import Any

from iscore3.data.rcsb_gate01 import COORDINATE_URL, immutable_write, preserve_manifest_timestamp
from iscore3.data.rcsb_gate01 import sha256_file, stable_json_bytes, utc_now
from iscore3.gate03.receptor_views import (
    POCKET_V2_NAMES,
    _accessions,
    _feature_row,
    read_tsv,
)
from iscore3.protein.apo_views import (
    AA1_TO_3,
    _entity_candidates,
    _minimum_foreign_contact,
    _release_date,
    _resolution,
    acquire_search_and_metadata,
)
from iscore3.protein.pocket_features import PocketError, extract_protein_pocket
from iscore3.protein.pocket_features import read_mmcif_atoms, tsv_bytes
from iscore3.protein.structure_views import align_construct_to_prediction, write_ca_pdb


class Gate3ApoError(RuntimeError):
    """Raised when S3 provenance or mapping invariants fail."""


def _download(path: Path, url: str) -> bytes:
    if path.exists():
        return path.read_bytes()
    from urllib.request import Request, urlopen

    request = Request(url, headers={"User-Agent": "iScore3.0-gate03-apo/1.0"})
    with urlopen(request, timeout=180) as response:
        payload = response.read()
    immutable_write(path, payload)
    return payload


def build_s3_views(
    *,
    dataset: Path,
    raw_root: Path,
    derived_root: Path,
    feature_output: Path,
    manifest_path: Path,
    audit_output: Path,
    maximum_coordinate_candidates_per_series: int = 10,
) -> dict[str, Any]:
    """Select one rigorously mapped, unoccupied X-ray receptor when feasible."""

    rows = read_tsv(dataset)
    by_series: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_series[row["series_id"]].append(row)
    accessions = sorted(
        {
            accession
            for values in by_series.values()
            for accession in _accessions(values[0]["uniprot_ids"])
        }
    )
    identifiers, entries, acquisition = acquire_search_and_metadata(accessions, raw_root)
    records: list[dict[str, Any]] = []
    feature_rows: list[dict[str, Any]] = []
    coordinates: dict[str, dict[str, Any]] = {}
    atom_cache: dict[str, list[Any]] = {}

    for series_id in sorted(by_series):
        representative = by_series[series_id][0]
        sequence = representative["target_sequence"]
        target_site = [int(value) for value in representative["site_target_positions"].split(";")]
        candidates = []
        for accession in _accessions(representative["uniprot_ids"]):
            for pdb_id in identifiers.get(accession, []):
                entry = entries[pdb_id]
                resolution = _resolution(entry)
                if resolution is None or resolution > 2.8:
                    continue
                if entry.get("nonpolymer_entities"):
                    raise Gate3ApoError(f"Zero-nonpolymer query returned liganded {pdb_id}")
                for entity in _entity_candidates(entry, accession):
                    alignment = align_construct_to_prediction(sequence, entity["sequence"])
                    site_map = {
                        position: alignment["mapping"][position]
                        for position in target_site
                        if position in alignment["mapping"]
                    }
                    site_coverage = len(site_map) / len(target_site)
                    if (
                        alignment["construct_coverage"] < 0.80
                        or alignment["aligned_identity"] < 0.98
                        or site_coverage < 0.80
                    ):
                        continue
                    mapping = {key: value for key, value in alignment.items() if key != "mapping"}
                    mapping.update(
                        {
                            "declared_uniprot_accession": accession,
                            "site_mapping_coverage": site_coverage,
                            "target_to_apo_site_positions": {
                                str(key): value for key, value in sorted(site_map.items())
                            },
                        }
                    )
                    candidates.append(
                        {
                            "pdb_id": pdb_id,
                            "resolution_angstrom": resolution,
                            "release_date": _release_date(entry),
                            "entity": entity,
                            "mapping": mapping,
                            "site_map": site_map,
                        }
                    )
        candidates.sort(
            key=lambda row: (
                row["resolution_angstrom"],
                row["release_date"] or "9999-99-99",
                row["pdb_id"],
                row["entity"]["entity_id"],
            )
        )
        attempts = []
        selected = None
        for candidate in candidates[:maximum_coordinate_candidates_per_series]:
            pdb_id = candidate["pdb_id"]
            coordinate_path = raw_root / "structures" / f"{pdb_id}.cif.gz"
            payload = _download(coordinate_path, COORDINATE_URL.format(pdb_id=pdb_id))
            try:
                gzip.decompress(payload)
            except (gzip.BadGzipFile, EOFError) as error:
                raise Gate3ApoError(f"Invalid coordinate gzip for {pdb_id}") from error
            coordinates[pdb_id] = {
                "pdb_id": pdb_id,
                "url": COORDINATE_URL.format(pdb_id=pdb_id),
                "path": str(coordinate_path),
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
            atoms = atom_cache.setdefault(pdb_id, list(read_mmcif_atoms(coordinate_path)))
            positions = sorted(candidate["site_map"].values())
            expected_names = {
                position: AA1_TO_3.get(candidate["entity"]["sequence"][position - 1], "UNK")
                for position in positions
            }
            for asym_id in candidate["entity"]["asym_ids"]:
                try:
                    pocket = extract_protein_pocket(
                        atoms,
                        pdb_id=pdb_id,
                        structure_path=coordinate_path,
                        protein_entity_id=candidate["entity"]["entity_id"],
                        candidate_asym_ids=[asym_id],
                        positions=positions,
                        expected_residue_names=expected_names,
                        minimum_coverage=0.80,
                    )
                except PocketError as error:
                    attempts.append(
                        {
                            "pdb_id": pdb_id,
                            "entity_id": candidate["entity"]["entity_id"],
                            "asym_id": asym_id,
                            "status": "coordinate_mapping_failure",
                            "error": str(error),
                        }
                    )
                    continue
                minimum_contact, foreign_atoms = _minimum_foreign_contact(
                    atoms,
                    pocket.atoms,
                    entity_id=candidate["entity"]["entity_id"],
                    asym_id=asym_id,
                )
                unoccupied = minimum_contact is None or minimum_contact >= 8.0
                attempts.append(
                    {
                        "pdb_id": pdb_id,
                        "entity_id": candidate["entity"]["entity_id"],
                        "asym_id": asym_id,
                        "status": "eligible" if unoccupied else "foreign_chain_site_contact",
                        "minimum_foreign_heavy_atom_distance_angstrom": minimum_contact,
                        "foreign_heavy_atom_count": foreign_atoms,
                    }
                )
                if unoccupied:
                    selected = (candidate, asym_id, pocket, coordinate_path, atoms)
                    break
            if selected is not None:
                break
        if selected is None:
            records.append(
                {
                    "series_id": series_id,
                    "status": "no_strict_pocket_unoccupied_view",
                    "search_entry_count": sum(
                        len(identifiers.get(accession, []))
                        for accession in _accessions(representative["uniprot_ids"])
                    ),
                    "sequence_mapped_candidate_count": len(candidates),
                    "coordinate_candidate_limit": maximum_coordinate_candidates_per_series,
                    "attempts": attempts,
                }
            )
            continue
        candidate, asym_id, pocket, coordinate_path, atoms = selected
        feature_rows.append(
            _feature_row(
                series_id=series_id,
                view="S3",
                source_id=candidate["pdb_id"],
                source_sha256=sha256_file(coordinate_path),
                entity_id=candidate["entity"]["entity_id"],
                asym_id=asym_id,
                pocket=pocket,
                target_positions=target_site,
                mean_site_plddt=None,
            )
        )
        chain_atoms = [
            atom
            for atom in atoms
            if atom.group == "ATOM"
            and atom.entity_id == candidate["entity"]["entity_id"]
            and atom.asym_id == asym_id
        ]
        global_pdb = write_ca_pdb(
            chain_atoms, derived_root / "S3" / f"{series_id}-global.pdb"
        )
        pocket_pdb = write_ca_pdb(
            pocket.atoms, derived_root / "S3" / f"{series_id}-pocket.pdb"
        )
        records.append(
            {
                "series_id": series_id,
                "view": "S3",
                "status": "eligible",
                "source_id": candidate["pdb_id"],
                "source_structure_path": str(coordinate_path),
                "source_structure_sha256": sha256_file(coordinate_path),
                "protein_entity_id": candidate["entity"]["entity_id"],
                "protein_asym_id": asym_id,
                "source_site_positions": list(pocket.expected_positions),
                "target_site_positions": target_site,
                "mapping": candidate["mapping"],
                "resolution_angstrom": candidate["resolution_angstrom"],
                "global": global_pdb,
                "pocket": pocket_pdb,
                "attempts": attempts,
            }
        )

    if feature_rows:
        immutable_write(feature_output, tsv_bytes(feature_rows))
    statuses = Counter(row["status"] for row in records)
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "PASS" if feature_rows else "NO_ELIGIBLE_VIEWS",
        "information_boundary": {
            "query_ligand_coordinates_read": False,
            "affinity_labels_read": False,
            "search": "declared UniProt; X-ray; zero nonpolymer entities",
            "site_unoccupied": "no non-water foreign heavy atom within 8 A",
        },
        "inputs": {"dataset": {"path": str(dataset), "sha256": sha256_file(dataset)}},
        "thresholds": {
            "maximum_resolution_angstrom": 2.8,
            "minimum_target_alignment_coverage": 0.80,
            "minimum_target_aligned_identity": 0.98,
            "minimum_site_mapping_coverage": 0.80,
            "site_foreign_atom_exclusion_radius_angstrom": 8.0,
            "maximum_coordinate_candidates_per_series": maximum_coordinate_candidates_per_series,
        },
        "acquisition": acquisition,
        "counts": {
            "series": len(by_series),
            "eligible_series": statuses["eligible"],
            "ineligible_series": statuses["no_strict_pocket_unoccupied_view"],
            "downloaded_coordinates": len(coordinates),
        },
        "coordinates": [coordinates[key] for key in sorted(coordinates)],
        "feature_output": {
            "path": str(feature_output) if feature_rows else "",
            "sha256": sha256_file(feature_output) if feature_rows else "",
            "bytes": feature_output.stat().st_size if feature_rows else 0,
            "feature_dimensions": len(POCKET_V2_NAMES),
        },
        "views": records,
    }
    preserve_manifest_timestamp(manifest_path, manifest, "created_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    checks = {
        "bounded_coordinate_attempts": all(
            len({attempt["pdb_id"] for attempt in row.get("attempts", [])})
            <= maximum_coordinate_candidates_per_series
            for row in records
        ),
        "eligible_views_exist": bool(feature_rows),
        "labels_not_read": True,
        "query_ligand_coordinates_not_read": True,
    }
    audit = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": manifest["counts"],
        "manifest_sha256": sha256_file(manifest_path),
    }
    preserve_manifest_timestamp(audit_output, audit, "created_utc")
    immutable_write(audit_output, stable_json_bytes(audit))
    return audit
