"""Strict pocket-unoccupied experimental receptor sensitivity views."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np

from iscore3.artifacts import (
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
    utc_now,
)
from iscore3.protein.rcsb_client import (
    COORDINATE_URL,
    DATA_URL,
    ENTRY_QUERY,
    SEARCH_URL,
    chunks,
    post_json,
)
from iscore3.protein.pocket_features import (
    FEATURE_NAMES,
    PocketError,
    extract_protein_pocket,
    pocket_feature_dict,
    read_mmcif_atoms,
    read_tsv,
    tsv_bytes,
)
from iscore3.protein.structure_views import AA1_TO_3, align_construct_to_prediction


USER_AGENT = "iScore3.0-receptor-apo/1.0 (scientific provenance audit)"
WATER_NAMES = {"HOH", "DOD", "WAT"}


class ApoViewError(RuntimeError):
    """Raised when the S3 acquisition or mapping contract is violated."""


def _download(path: Path, url: str) -> bytes:
    if path.exists():
        return path.read_bytes()
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=180) as response:
        payload = response.read()
    immutable_write(path, payload)
    return payload


def _search_payload(accession: str) -> dict[str, Any]:
    nodes = []
    for attribute, operator, value in (
        (
            "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            "exact_match",
            accession,
        ),
        (
            "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
            "exact_match",
            "UniProt",
        ),
        ("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        ("rcsb_entry_info.nonpolymer_entity_count", "equals", 0),
    ):
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": attribute,
                    "operator": operator,
                    "value": value,
                },
            }
        )
    return {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {"return_all_hits": True, "results_verbosity": "compact"},
    }


def _post_search(payload: Mapping[str, Any], attempts: int = 5) -> dict[str, Any]:
    body = json.dumps(payload, separators=(",", ":"), allow_nan=False).encode("utf-8")
    last_error: Exception | None = None
    for attempt in range(attempts):
        request = Request(
            SEARCH_URL,
            data=body,
            headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
            method="POST",
        )
        try:
            with urlopen(request, timeout=120) as response:
                raw = response.read()
            if not raw:
                return {"total_count": 0, "result_set": []}
            value = json.loads(raw.decode("utf-8"))
            if not isinstance(value, dict):
                raise ApoViewError("Non-object RCSB search response")
            return value
        except HTTPError as error:
            if error.code == 204:
                return {"total_count": 0, "result_set": []}
            last_error = error
        except (
            URLError,
            TimeoutError,
            ConnectionResetError,
            json.JSONDecodeError,
        ) as error:
            last_error = error
        if attempt + 1 < attempts:
            time.sleep(2**attempt)
    raise ApoViewError(
        f"RCSB apo search failed after {attempts} attempts: {last_error}"
    )


def acquire_search_and_metadata(
    accessions: Sequence[str], raw_root: Path, *, batch_size: int = 50
) -> tuple[dict[str, list[str]], dict[str, Mapping[str, Any]], dict[str, Any]]:
    search_root = raw_root / "search"
    identifiers_by_accession: dict[str, list[str]] = {}
    search_records = []
    all_identifiers: set[str] = set()
    for accession in sorted(set(accessions)):
        path = search_root / f"{accession}.json"
        payload = _search_payload(accession)
        if path.exists():
            response = json.loads(path.read_text(encoding="utf-8"))
        else:
            response = _post_search(payload)
            immutable_write(path, stable_json_bytes(response))
        identifiers = sorted(
            {
                str(item if isinstance(item, str) else item["identifier"]).upper()
                for item in response.get("result_set", [])
            }
        )
        total = int(response.get("total_count", len(identifiers)))
        if total != len(identifiers):
            raise ApoViewError(f"Incomplete RCSB apo search for {accession}: {total}")
        identifiers_by_accession[accession] = identifiers
        all_identifiers.update(identifiers)
        search_records.append(
            {
                "uniprot_accession": accession,
                "request": payload,
                "response_path": str(path),
                "response_sha256": sha256_file(path),
                "response_bytes": path.stat().st_size,
                "entry_count": len(identifiers),
            }
        )

    entries: dict[str, Mapping[str, Any]] = {}
    metadata_records = []
    for index, batch in enumerate(chunks(sorted(all_identifiers), batch_size)):
        path = raw_root / "metadata" / f"entries-{index:04d}.json"
        if path.exists():
            response = json.loads(path.read_text(encoding="utf-8"))
        else:
            response = post_json(
                DATA_URL, {"query": ENTRY_QUERY, "variables": {"ids": batch}}
            )
            immutable_write(path, stable_json_bytes(response))
        returned = response.get("data", {}).get("entries") or []
        for entry in returned:
            if entry:
                entries[str(entry["rcsb_id"]).upper()] = entry
        metadata_records.append(
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "requested": len(batch),
                "returned": len(returned),
            }
        )
    missing = sorted(all_identifiers.difference(entries))
    if missing:
        raise ApoViewError(f"RCSB Data API omitted apo candidates: {missing[:20]}")
    provenance = {
        "search_url": SEARCH_URL,
        "data_url": DATA_URL,
        "snapshot_semantics": "mutable_service_responses_cached_immutably",
        "query_policy": "same exact UniProt; X-ray; zero nonpolymer entities",
        "search_records": search_records,
        "metadata_records": metadata_records,
        "unique_candidate_entries": len(all_identifiers),
    }
    return identifiers_by_accession, entries, provenance


def _resolution(entry: Mapping[str, Any]) -> float | None:
    values = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
    parsed = [float(value) for value in values if value is not None]
    return min(parsed) if parsed else None


def _release_date(entry: Mapping[str, Any]) -> str:
    return str(
        (entry.get("rcsb_accession_info") or {}).get("initial_release_date") or ""
    )


def _entity_candidates(
    entry: Mapping[str, Any], accession: str
) -> list[dict[str, Any]]:
    result = []
    for entity in entry.get("polymer_entities") or []:
        polymer = entity.get("entity_poly") or {}
        if polymer.get("rcsb_entity_polymer_type") != "Protein":
            continue
        identifiers = entity.get("rcsb_polymer_entity_container_identifiers") or {}
        accessions = {
            str(reference.get("database_accession") or "").upper()
            for reference in identifiers.get("reference_sequence_identifiers") or []
            if str(reference.get("database_name") or "").upper() == "UNIPROT"
        }
        if accession.upper() not in accessions:
            continue
        sequence = "".join(
            str(polymer.get("pdbx_seq_one_letter_code_can") or "").split()
        ).upper()
        if not sequence or not sequence.isalpha():
            continue
        result.append(
            {
                "entity_id": str(identifiers.get("entity_id") or ""),
                "asym_ids": [str(value) for value in identifiers.get("asym_ids") or []],
                "sequence": sequence,
            }
        )
    return result


def _minimum_foreign_contact(
    all_atoms: Sequence[Any],
    pocket_atoms: Sequence[Any],
    *,
    entity_id: str,
    asym_id: str,
) -> tuple[float | None, int]:
    pocket_xyz = np.asarray([atom.xyz for atom in pocket_atoms], dtype=np.float64)
    foreign = [
        atom
        for atom in all_atoms
        if not (atom.entity_id == entity_id and atom.asym_id == asym_id)
        and atom.residue_name.upper() not in WATER_NAMES
        and atom.element.upper() not in {"H", "D"}
        and atom.group in {"ATOM", "HETATM"}
    ]
    if not foreign:
        return None, 0
    minimum = math.inf
    for start in range(0, len(foreign), 4096):
        coordinates = np.asarray(
            [atom.xyz for atom in foreign[start : start + 4096]], dtype=np.float64
        )
        distance = np.sqrt(
            np.sum((coordinates[:, None, :] - pocket_xyz[None, :, :]) ** 2, axis=2)
        )
        minimum = min(minimum, float(np.min(distance)))
    return minimum, len(foreign)


def _s3_row(
    template: Mapping[str, str],
    *,
    pdb_id: str,
    structure_sha256: str,
    entity_id: str,
    asym_id: str,
    positions: Sequence[int],
    present: Sequence[int],
    missing: Sequence[int],
    features: Mapping[str, float],
) -> dict[str, Any]:
    row: dict[str, Any] = dict(template)
    row.update(
        {
            "view_id": f"{template['observation_id']}:S3",
            "mapping_tier": "S3",
            "feature_structure_pdb_id": pdb_id,
            "feature_structure_sha256": structure_sha256,
            "feature_protein_entity_id": entity_id,
            "feature_protein_asym_id": asym_id,
            "pocket_positions_label_seq_id": ";".join(map(str, positions)),
            "present_positions_label_seq_id": ";".join(map(str, present)),
            "missing_positions_label_seq_id": ";".join(map(str, missing)),
            "receptor_conformation_source": "rcsb_pocket_unoccupied_experimental_fixed_per_construct",
            "query_holo_receptor_privilege": False,
            "historical_reference_ligand_used_only_for_site": True,
            "query_ligand_coordinates_read": False,
        }
    )
    row.update(features)
    return row


def build_apo_views(
    *,
    pilot: Path,
    strict_pockets: Path,
    sites: Path,
    raw_root: Path,
    output: Path,
    manifest_path: Path,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    all_rows = read_tsv(pilot)
    supervised = [row for row in all_rows if row["role"] == "supervised_s0"]
    by_group: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in supervised:
        by_group[row["construct_group_id"]].append(row)
    definitions = {
        definition["construct_group_id"]: definition
        for definition in json.loads(sites.read_text(encoding="utf-8"))["definitions"]
    }
    templates = {
        row["observation_id"]: row
        for row in read_tsv(strict_pockets)
        if row["mapping_tier"] == "S1"
    }
    accessions = sorted({rows[0]["uniprot_accession"] for rows in by_group.values()})
    identifiers_by_accession, entries, acquisition = acquire_search_and_metadata(
        accessions, raw_root
    )
    policy = config["apo_view"]
    records = []
    output_rows = []
    coordinate_records: dict[str, dict[str, Any]] = {}
    atom_cache: dict[str, list[Any]] = {}

    for group_id in sorted(by_group):
        representative = by_group[group_id][0]
        accession = representative["uniprot_accession"]
        site_positions = [
            int(value) for value in definitions[group_id]["positions_label_seq_id"]
        ]
        candidates = []
        for pdb_id in identifiers_by_accession.get(accession, []):
            entry = entries[pdb_id]
            resolution = _resolution(entry)
            if resolution is None or resolution > 2.8:
                continue
            if entry.get("nonpolymer_entities"):
                raise ApoViewError(
                    f"Zero-nonpolymer search returned nonpolymer entry {pdb_id}"
                )
            for entity in _entity_candidates(entry, accession):
                alignment = align_construct_to_prediction(
                    representative["construct_sequence"], entity["sequence"]
                )
                site_map = {
                    position: alignment["mapping"][position]
                    for position in site_positions
                    if position in alignment["mapping"]
                }
                site_coverage = len(site_map) / len(site_positions)
                mapping = {
                    key: value for key, value in alignment.items() if key != "mapping"
                }
                mapping.update(
                    {
                        "site_mapping_coverage": site_coverage,
                        "mapped_site_positions": {
                            str(key): value for key, value in sorted(site_map.items())
                        },
                    }
                )
                if (
                    alignment["construct_coverage"]
                    < float(policy["minimum_construct_alignment_coverage"])
                    or alignment["aligned_identity"]
                    < float(policy["minimum_construct_aligned_identity"])
                    or site_coverage < float(policy["minimum_site_mapping_coverage"])
                ):
                    continue
                candidates.append(
                    {
                        "pdb_id": pdb_id,
                        "entry": entry,
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
        for candidate in candidates:
            pdb_id = candidate["pdb_id"]
            coordinate_path = raw_root / "structures" / f"{pdb_id}.cif.gz"
            payload = _download(coordinate_path, COORDINATE_URL.format(pdb_id=pdb_id))
            try:
                gzip.decompress(payload)
            except (gzip.BadGzipFile, EOFError) as error:
                raise ApoViewError(f"Invalid coordinate gzip for {pdb_id}") from error
            coordinate_records[pdb_id] = {
                "pdb_id": pdb_id,
                "path": str(coordinate_path),
                "url": COORDINATE_URL.format(pdb_id=pdb_id),
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
            atoms = atom_cache.setdefault(pdb_id, read_mmcif_atoms(coordinate_path))
            target_positions = sorted(candidate["site_map"].values())
            expected_names = {
                position: AA1_TO_3.get(
                    candidate["entity"]["sequence"][position - 1], "UNK"
                )
                for position in target_positions
            }
            for asym_id in candidate["entity"]["asym_ids"]:
                try:
                    pocket = extract_protein_pocket(
                        atoms,
                        pdb_id=pdb_id,
                        structure_path=coordinate_path,
                        protein_entity_id=candidate["entity"]["entity_id"],
                        candidate_asym_ids=[asym_id],
                        positions=target_positions,
                        expected_residue_names=expected_names,
                        minimum_coverage=float(policy["minimum_site_mapping_coverage"]),
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
                unoccupied = minimum_contact is None or minimum_contact >= float(
                    policy["site_nonpolymer_exclusion_radius_angstrom"]
                )
                attempts.append(
                    {
                        "pdb_id": pdb_id,
                        "entity_id": candidate["entity"]["entity_id"],
                        "asym_id": asym_id,
                        "status": (
                            "eligible" if unoccupied else "foreign_chain_site_contact"
                        ),
                        "resolution_angstrom": candidate["resolution_angstrom"],
                        "mapping": candidate["mapping"],
                        "foreign_heavy_atom_count": foreign_atoms,
                        "minimum_foreign_heavy_atom_distance_angstrom": minimum_contact,
                    }
                )
                if not unoccupied:
                    continue
                selected = (candidate, asym_id, pocket, coordinate_path)
                break
            if selected is not None:
                break
        if selected is None:
            records.append(
                {
                    "construct_group_id": group_id,
                    "uniprot_accession": accession,
                    "status": "no_strict_pocket_unoccupied_view",
                    "search_entry_count": len(
                        identifiers_by_accession.get(accession, [])
                    ),
                    "sequence_mapped_candidate_count": len(candidates),
                    "attempts": attempts,
                }
            )
            continue
        candidate, asym_id, pocket, coordinate_path = selected
        features = pocket_feature_dict(pocket)
        for row in by_group[group_id]:
            output_rows.append(
                _s3_row(
                    templates[row["observation_id"]],
                    pdb_id=candidate["pdb_id"],
                    structure_sha256=sha256_file(coordinate_path),
                    entity_id=candidate["entity"]["entity_id"],
                    asym_id=asym_id,
                    positions=pocket.expected_positions,
                    present=pocket.present_positions,
                    missing=pocket.missing_positions,
                    features=features,
                )
            )
        records.append(
            {
                "construct_group_id": group_id,
                "construct_sha256": representative["construct_sha256"],
                "uniprot_accession": accession,
                "status": "eligible",
                "source_pdb_id": candidate["pdb_id"],
                "source_structure_path": str(coordinate_path),
                "source_structure_sha256": sha256_file(coordinate_path),
                "protein_entity_id": candidate["entity"]["entity_id"],
                "protein_asym_id": asym_id,
                "resolution_angstrom": candidate["resolution_angstrom"],
                "mapping": candidate["mapping"],
                "attempts": attempts,
            }
        )

    if output_rows:
        immutable_write(output, tsv_bytes(output_rows))
    statuses = Counter(record["status"] for record in records)
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "PASS" if output_rows else "NO_ELIGIBLE_VIEWS",
        "information_boundary": {
            "query_ligand_coordinates_read": False,
            "site_source": "historical-reference residue positions only",
            "candidate_search": "same exact UniProt, X-ray, zero nonpolymer entities",
            "pocket_unoccupied_definition": (
                "no non-water foreign polymer/nonpolymer heavy atom within frozen 8-A site radius"
            ),
        },
        "inputs": {
            "pilot": {"path": str(pilot), "sha256": sha256_file(pilot)},
            "strict_pockets": {
                "path": str(strict_pockets),
                "sha256": sha256_file(strict_pockets),
            },
            "sites": {"path": str(sites), "sha256": sha256_file(sites)},
        },
        "thresholds": policy,
        "acquisition": acquisition,
        "counts": {
            "construct_groups": len(by_group),
            "eligible_construct_groups": statuses["eligible"],
            "ineligible_construct_groups": statuses["no_strict_pocket_unoccupied_view"],
            "S3_observation_rows": len(output_rows),
            "downloaded_coordinate_entries": len(coordinate_records),
        },
        "coordinates": [coordinate_records[key] for key in sorted(coordinate_records)],
        "output": {
            "path": str(output) if output_rows else "",
            "sha256": sha256_file(output) if output_rows else "",
            "bytes": output.stat().st_size if output_rows else 0,
            "feature_schema": list(FEATURE_NAMES),
        },
        "views": records,
    }
    preserve_manifest_timestamp(manifest_path, manifest, "created_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    return manifest
