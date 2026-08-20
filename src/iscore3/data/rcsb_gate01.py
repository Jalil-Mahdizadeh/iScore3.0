"""RCSB-derived, entry-linked Kd pilot acquisition for Gate-0/1.

This module treats the RCSB Search/Data APIs as mutable services. Every response
is cached immutably and SHA-256-addressed. The selected pilot is deliberately
small: exact structure-linked Kd annotations, one protein entity, an explicit
co-complex ligand component, and exact construct-sequence grouping.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
import csv
import gzip
import hashlib
from http.client import RemoteDisconnected
import json
import math
from pathlib import Path
import time
from typing import Any, Iterable, Iterator, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/graphql"
COORDINATE_URL = "https://files.rcsb.org/download/{pdb_id}.cif.gz"
USER_AGENT = "iScore3-Gate01/0.1 (scientific-reproducibility; RCSB-cached-requests)"

ENTRY_QUERY = r"""
query Gate01Entries($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_accession_info { deposit_date initial_release_date revision_date }
    rcsb_binding_affinity { comp_id type value unit }
    exptl { method }
    rcsb_entry_info { resolution_combined }
    rcsb_primary_citation { title pdbx_database_id_DOI pdbx_database_id_PubMed }
    polymer_entities {
      rcsb_id
      entity_poly { pdbx_seq_one_letter_code_can rcsb_entity_polymer_type }
      rcsb_polymer_entity { pdbx_description pdbx_mutation pdbx_fragment }
      rcsb_polymer_entity_container_identifiers {
        entity_id asym_ids auth_asym_ids
        reference_sequence_identifiers { database_accession database_name }
      }
      rcsb_entity_source_organism { ncbi_scientific_name ncbi_taxonomy_id }
    }
    nonpolymer_entities {
      rcsb_id
      rcsb_nonpolymer_entity_container_identifiers {
        entity_id asym_ids auth_asym_ids nonpolymer_comp_id
      }
      nonpolymer_comp {
        chem_comp { id type formula_weight formula }
        rcsb_chem_comp_descriptor { SMILES SMILES_stereo InChI InChIKey }
      }
    }
  }
}
""".strip()

ALLOWED_ELEMENTS = {"C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "H", "B", "Si"}


class RcsbPilotError(RuntimeError):
    """Raised when an acquisition or curation invariant fails."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def stable_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")


def immutable_write(path: Path, payload: bytes) -> str:
    """Write once, or verify that an existing snapshot is byte-identical."""

    path.parent.mkdir(parents=True, exist_ok=True)
    expected = sha256_bytes(payload)
    if path.exists():
        observed = sha256_file(path)
        if observed != expected:
            raise RcsbPilotError(f"Refusing to replace nonidentical immutable file: {path}")

        return observed
    temporary = path.with_name(f".{path.name}.partial")
    if temporary.exists():
        raise RcsbPilotError(f"Stale partial output requires review: {temporary}")
    with temporary.open("xb") as handle:
        handle.write(payload)
        handle.flush()
    temporary.replace(path)
    return expected
def preserve_manifest_timestamp(path: Path, manifest: dict[str, Any], field: str) -> None:
    """Keep write-once manifests byte-stable when a command is rerun."""

    if not path.exists():
        return
    previous = json.loads(path.read_text(encoding="utf-8"))
    if field not in previous:
        raise RcsbPilotError(f"Existing manifest lacks required timestamp {field}: {path}")
    manifest[field] = previous[field]



def post_json(url: str, payload: Mapping[str, Any], *, attempts: int = 4) -> dict[str, Any]:
    body = json.dumps(payload, separators=(",", ":"), allow_nan=False).encode("utf-8")
    request = Request(
        url,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        method="POST",
    )
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            with urlopen(request, timeout=120) as response:
                result = json.loads(response.read().decode("utf-8"))
            if not isinstance(result, dict):
                raise RcsbPilotError(f"Non-object JSON response from {url}")
            if result.get("errors"):
                raise RcsbPilotError(f"Graph/API errors from {url}: {result['errors']}")
            return result
        except (
            HTTPError,
            URLError,
            TimeoutError,
            RemoteDisconnected,
            ConnectionResetError,
            json.JSONDecodeError,
        ) as error:
            last_error = error
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RcsbPilotError(f"Request failed after {attempts} attempts: {last_error}")


def chunks(values: Sequence[str], size: int) -> Iterator[list[str]]:
    if size <= 0:
        raise ValueError("Chunk size must be positive")
    for start in range(0, len(values), size):
        yield list(values[start : start + size])


def acquire_metadata(
    raw_root: Path,
    manifest_path: Path,
    *,
    endpoint: str = "Kd",
    batch_size: int = 100,
) -> dict[str, Any]:
    """Cache the complete current RCSB entry set for one affinity endpoint."""

    raw_root = raw_root.resolve()
    search_path = raw_root / f"search-{endpoint.lower()}.json"
    search_payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_binding_affinity.type",
                "operator": "exact_match",
                "value": endpoint,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 10000}},
    }
    if search_path.exists():
        search_response = json.loads(search_path.read_text(encoding="utf-8"))
    else:
        search_response = post_json(SEARCH_URL, search_payload)
        immutable_write(search_path, stable_json_bytes(search_response))
    identifiers = sorted(
        {str(item["identifier"]).upper() for item in search_response.get("result_set", [])}
    )
    total_count = int(search_response.get("total_count", -1))
    if not identifiers or total_count != len(identifiers):
        raise RcsbPilotError(
            f"RCSB search pagination incomplete: total={total_count}, identifiers={len(identifiers)}"
        )

    files: list[dict[str, Any]] = [
        {
            "path": str(search_path),
            "sha256": sha256_file(search_path),
            "bytes": search_path.stat().st_size,
            "kind": "search_response",
        }
    ]
    returned_ids: set[str] = set()
    for index, batch in enumerate(chunks(identifiers, batch_size)):
        path = raw_root / "metadata" / f"entries-{index:04d}.json"
        if path.exists():
            response = json.loads(path.read_text(encoding="utf-8"))
        else:
            response = post_json(DATA_URL, {"query": ENTRY_QUERY, "variables": {"ids": batch}})
            immutable_write(path, stable_json_bytes(response))
        entries = response.get("data", {}).get("entries") or []
        returned_ids.update(str(entry["rcsb_id"]).upper() for entry in entries if entry)
        files.append(
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "kind": "metadata_response",
                "requested_ids": len(batch),
                "returned_entries": len(entries),
            }
        )
    missing = sorted(set(identifiers) - returned_ids)
    if missing:
        raise RcsbPilotError(f"RCSB Data API omitted {len(missing)} searched entries: {missing[:20]}")

    manifest = {
        "schema_version": 1,
        "source": "RCSB PDB Search API and Data API",
        "source_urls": [SEARCH_URL, DATA_URL],
        "snapshot_semantics": "mutable_service_response_cached_immutably",
        "acquired_utc": utc_now(),
        "endpoint": endpoint,
        "entry_count": len(identifiers),
        "batch_size": batch_size,
        "search_request": search_payload,
        "graphql_query_sha256": sha256_bytes(ENTRY_QUERY.encode("utf-8")),
        "files": files,
    }
    preserve_manifest_timestamp(manifest_path.resolve(), manifest, "acquired_utc")
    immutable_write(manifest_path.resolve(), stable_json_bytes(manifest))
    return manifest


def iter_cached_entries(raw_root: Path) -> Iterator[dict[str, Any]]:
    paths = sorted((raw_root / "metadata").glob("entries-*.json"))
    if not paths:
        raise RcsbPilotError(f"No cached RCSB metadata batches under {raw_root}")
    for path in paths:
        response = json.loads(path.read_text(encoding="utf-8"))
        for entry in response.get("data", {}).get("entries") or []:
            if entry:
                yield entry


def _first(values: Any) -> Any | None:
    return values[0] if isinstance(values, list) and values else None


def _date(value: Any) -> str:
    return str(value or "").split("T", 1)[0]


def _sequence(value: Any) -> str:
    return "".join(str(value or "").split()).replace("(", "").replace(")", "").upper()


def _reference_accessions(identifiers: Mapping[str, Any]) -> list[str]:
    return sorted(
        {
            str(item["database_accession"])
            for item in identifiers.get("reference_sequence_identifiers") or []
            if item and str(item.get("database_name", "")).lower() == "uniprot"
        }
    )


@dataclass(frozen=True, slots=True)
class Candidate:
    pdb_id: str
    ligand_comp_id: str
    endpoint: str
    value_nm: float
    pKd: float
    canonical_smiles: str
    inchikey: str
    ligand_heavy_atoms: int
    ligand_formula_weight: float
    protein_entity_id: str
    protein_asym_ids: str
    protein_auth_asym_ids: str
    construct_sequence: str
    construct_sha256: str
    construct_length: int
    uniprot_accession: str
    protein_description: str
    mutation_annotation: str
    organism: str
    taxonomy_id: str
    ligand_entity_id: str
    ligand_asym_ids: str
    ligand_auth_asym_ids: str
    resolution_angstrom: float
    experimental_method: str
    deposit_date: str
    release_date: str
    revision_date: str
    citation_doi: str
    citation_pubmed: str
    citation_title: str
    source_record_url: str


def candidate_from_entry(entry: Mapping[str, Any]) -> tuple[Candidate | None, str | None]:
    """Apply deliberately strict, label-blind candidate filters."""

    try:
        from rdkit import Chem
    except ImportError as error:  # pragma: no cover
        raise RcsbPilotError("RDKit is required for pilot ligand standardization") from error

    pdb_id = str(entry.get("rcsb_id", "")).upper()
    affinities = [
        item
        for item in entry.get("rcsb_binding_affinity") or []
        if item and str(item.get("type")) == "Kd"
    ]
    if len(affinities) != 1:
        return None, "not_exactly_one_kd_annotation"
    affinity = affinities[0]
    if str(affinity.get("unit")) != "nM":
        return None, "kd_unit_not_nm"
    try:
        value_nm = float(affinity["value"])
    except (KeyError, TypeError, ValueError):
        return None, "invalid_kd_value"
    if not math.isfinite(value_nm) or value_nm <= 0:
        return None, "nonpositive_or_nonfinite_kd"

    methods = sorted({str(item.get("method")) for item in entry.get("exptl") or [] if item})
    if methods != ["X-RAY DIFFRACTION"]:
        return None, "not_xray_only"
    resolutions = entry.get("rcsb_entry_info", {}).get("resolution_combined") or []
    if len(resolutions) != 1:
        return None, "missing_or_multiple_resolution"
    resolution = float(resolutions[0])
    if not math.isfinite(resolution) or resolution > 2.8:
        return None, "resolution_above_2p8"

    proteins = [
        entity
        for entity in entry.get("polymer_entities") or []
        if entity
        and str((entity.get("entity_poly") or {}).get("rcsb_entity_polymer_type")) == "Protein"
    ]
    if len(proteins) != 1:
        return None, "not_exactly_one_protein_entity"
    protein = proteins[0]
    entity_poly = protein.get("entity_poly") or {}
    sequence = _sequence(entity_poly.get("pdbx_seq_one_letter_code_can"))
    if not 50 <= len(sequence) <= 2000 or not set(sequence).issubset(set("ACDEFGHIKLMNPQRSTVWYUX")):
        return None, "invalid_construct_sequence"
    identifiers = protein.get("rcsb_polymer_entity_container_identifiers") or {}
    uniprot = _reference_accessions(identifiers)
    if len(uniprot) != 1:
        return None, "not_exactly_one_uniprot_accession"

    comp_id = str(affinity.get("comp_id", "")).upper()
    ligands = [
        entity
        for entity in entry.get("nonpolymer_entities") or []
        if str(
            (entity.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get(
                "nonpolymer_comp_id", ""
            )
        ).upper()
        == comp_id
    ]
    if len(ligands) != 1:
        return None, "affinity_ligand_entity_not_unique"
    ligand = ligands[0]
    ligand_ids = ligand.get("rcsb_nonpolymer_entity_container_identifiers") or {}
    component = ligand.get("nonpolymer_comp") or {}
    chem_comp = component.get("chem_comp") or {}
    descriptor = component.get("rcsb_chem_comp_descriptor") or {}
    smiles = descriptor.get("SMILES_stereo") or descriptor.get("SMILES")
    molecule = Chem.MolFromSmiles(str(smiles or ""))
    if molecule is None:
        return None, "ligand_smiles_parse_failure"
    if len(Chem.GetMolFrags(molecule)) != 1:
        return None, "ligand_disconnected"
    elements = {atom.GetSymbol() for atom in molecule.GetAtoms()}
    if "C" not in elements or not elements.issubset(ALLOWED_ELEMENTS):
        return None, "ligand_element_policy"
    if any(atom.GetAtomMapNum() or atom.GetIsotope() for atom in molecule.GetAtoms()):
        return None, "ligand_annotation_policy"
    heavy_atoms = int(molecule.GetNumHeavyAtoms())
    if not 6 <= heavy_atoms <= 60:
        return None, "ligand_heavy_atom_range"
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    formula_weight = float(chem_comp.get("formula_weight") or float("nan"))
    if not math.isfinite(formula_weight) or not 80 <= formula_weight <= 900:
        return None, "ligand_formula_weight_range"
    inchikey = str(descriptor.get("InChIKey") or "")
    if not inchikey:
        return None, "missing_inchikey"

    citation = entry.get("rcsb_primary_citation") or {}
    doi = str(citation.get("pdbx_database_id_DOI") or "")
    pubmed = str(citation.get("pdbx_database_id_PubMed") or "")
    if not doi and not pubmed:
        return None, "missing_primary_citation_identifier"
    accession = entry.get("rcsb_accession_info") or {}
    protein_meta = protein.get("rcsb_polymer_entity") or {}
    organism_row = _first(protein.get("rcsb_entity_source_organism")) or {}

    return (
        Candidate(
            pdb_id=pdb_id,
            ligand_comp_id=comp_id,
            endpoint="Kd",
            value_nm=value_nm,
            pKd=9.0 - math.log10(value_nm),
            canonical_smiles=canonical,
            inchikey=inchikey,
            ligand_heavy_atoms=heavy_atoms,
            ligand_formula_weight=formula_weight,
            protein_entity_id=str(identifiers.get("entity_id")),
            protein_asym_ids=";".join(sorted(str(x) for x in identifiers.get("asym_ids") or [])),
            protein_auth_asym_ids=";".join(
                sorted(str(x) for x in identifiers.get("auth_asym_ids") or [])
            ),
            construct_sequence=sequence,
            construct_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            construct_length=len(sequence),
            uniprot_accession=uniprot[0],
            protein_description=str(protein_meta.get("pdbx_description") or ""),
            mutation_annotation=str(protein_meta.get("pdbx_mutation") or ""),
            organism=str(organism_row.get("ncbi_scientific_name") or ""),
            taxonomy_id=str(organism_row.get("ncbi_taxonomy_id") or ""),
            ligand_entity_id=str(ligand_ids.get("entity_id")),
            ligand_asym_ids=";".join(sorted(str(x) for x in ligand_ids.get("asym_ids") or [])),
            ligand_auth_asym_ids=";".join(
                sorted(str(x) for x in ligand_ids.get("auth_asym_ids") or [])
            ),
            resolution_angstrom=resolution,
            experimental_method=methods[0],
            deposit_date=_date(accession.get("deposit_date")),
            release_date=_date(accession.get("initial_release_date")),
            revision_date=_date(accession.get("revision_date")),
            citation_doi=doi,
            citation_pubmed=pubmed,
            citation_title=str(citation.get("title") or ""),
            source_record_url=f"https://www.rcsb.org/structure/{pdb_id}",
        ),
        None,
    )


def _candidate_sort_key(candidate: Candidate) -> tuple[Any, ...]:
    return (
        candidate.release_date or "9999-99-99",
        candidate.resolution_angstrom,
        candidate.pdb_id,
    )


def _stable_member_key(candidate: Candidate, seed: int) -> str:
    return hashlib.sha256(
        f"{seed}\x1f{candidate.construct_sha256}\x1f{candidate.inchikey}\x1f{candidate.pdb_id}".encode(
            "utf-8"
        )
    ).hexdigest()


def select_pilot(
    raw_root: Path,
    output_tsv: Path,
    manifest_path: Path,
    *,
    min_supervised_per_construct: int = 8,
    max_supervised_per_construct: int = 20,
    max_constructs: int = 12,
    replicate_tolerance_pkd: float = 0.30,
    selection_seed: int = 20260820,
) -> dict[str, Any]:
    """Build an S0 pilot plus one quarantined site reference per construct."""

    rejections: Counter[str] = Counter()
    candidates: list[Candidate] = []
    for entry in iter_cached_entries(raw_root):
        candidate, reason = candidate_from_entry(entry)
        if reason:
            rejections[reason] += 1
        else:
            assert candidate is not None
            candidates.append(candidate)

    identity_groups: dict[tuple[str, str], list[Candidate]] = defaultdict(list)
    for candidate in candidates:
        identity_groups[(candidate.construct_sha256, candidate.inchikey)].append(candidate)
    deduplicated: list[Candidate] = []
    replicate_records: dict[str, dict[str, Any]] = {}
    for (construct_hash, inchikey), rows in sorted(identity_groups.items()):
        values = [row.pKd for row in rows]
        key = f"{construct_hash}:{inchikey}"
        if len(rows) > 1:
            replicate_records[key] = {
                "pdb_ids": sorted(row.pdb_id for row in rows),
                "pKd_values": sorted(values),
                "range_pKd": max(values) - min(values),
            }
        if max(values) - min(values) > replicate_tolerance_pkd:
            rejections["discordant_structure_linked_replicates"] += len(rows)
            continue
        chosen = min(rows, key=lambda row: (row.resolution_angstrom, *_candidate_sort_key(row)))
        deduplicated.append(chosen)
        rejections["concordant_replicate_rows_collapsed"] += len(rows) - 1

    by_construct: dict[str, list[Candidate]] = defaultdict(list)
    for candidate in deduplicated:
        by_construct[candidate.construct_sha256].append(candidate)
    eligible = [
        rows
        for rows in by_construct.values()
        if len(rows) - 1 >= min_supervised_per_construct
    ]
    eligible.sort(
        key=lambda rows: (
            -min(len(rows) - 1, max_supervised_per_construct),
            rows[0].construct_sha256,
        )
    )
    selected_groups = eligible[:max_constructs]

    output_rows: list[dict[str, Any]] = []
    for group_index, rows in enumerate(selected_groups, start=1):
        reference = min(rows, key=_candidate_sort_key)
        members = [row for row in rows if row != reference]
        members.sort(key=lambda row: _stable_member_key(row, selection_seed))
        members = members[:max_supervised_per_construct]
        group_id = f"construct-{group_index:03d}-{reference.construct_sha256[:12]}"
        for role, row in [("site_reference_only", reference), *[("supervised_s0", x) for x in members]]:
            record = asdict(row)
            record.update(
                {
                    "observation_id": f"RCSB-{row.pdb_id}-{row.ligand_comp_id}-KD",
                    "construct_group_id": group_id,
                    "role": role,
                    "mapping_tier": "S0",
                    "site_reference_pdb_id": reference.pdb_id,
                    "site_reference_ligand_comp_id": reference.ligand_comp_id,
                    "site_definition_policy": "reference_ligand_6A_then_exact_sequence_position_transfer",
                    "query_ligand_coordinates_allowed": False,
                    "holo_receptor_privilege": True,
                }
            )
            record["label_quarantined"] = role == "site_reference_only"
            if role == "site_reference_only":
                record["value_nm"] = ""
                record["pKd"] = ""
                record["label_quarantine_reason"] = "historical_ligand_defines_site_only"
            else:
                record["label_quarantine_reason"] = ""
            output_rows.append(record)
    if not output_rows:
        raise RcsbPilotError("Strict selection produced no eligible construct groups")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(output_rows[0])
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(output_rows)
    payload = buffer.getvalue().encode("utf-8")
    immutable_write(output_tsv, payload)

    supervised = [row for row in output_rows if row["role"] == "supervised_s0"]
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "selection_stage": "before_model_fitting",
        "source_snapshot_root": str(raw_root.resolve()),
        "filters": {
            "endpoint": "Kd",
            "unit": "nM",
            "exactly_one_affinity_annotation": True,
            "experimental_method": "X-RAY DIFFRACTION",
            "maximum_resolution_angstrom": 2.8,
            "exactly_one_protein_entity": True,
            "exactly_one_uniprot_accession": True,
            "ligand_heavy_atoms": [6, 60],
            "ligand_formula_weight": [80, 900],
            "allowed_elements": sorted(ALLOWED_ELEMENTS),
            "replicate_tolerance_pKd": replicate_tolerance_pkd,
            "min_supervised_per_construct": min_supervised_per_construct,
            "max_supervised_per_construct": max_supervised_per_construct,
            "max_constructs": max_constructs,
            "selection_seed": selection_seed,
        },
        "counts": {
            "api_entries": sum(1 for _ in iter_cached_entries(raw_root)),
            "strict_candidates_before_deduplication": len(candidates),
            "deduplicated_candidates": len(deduplicated),
            "eligible_construct_groups": len(eligible),
            "selected_construct_groups": len(selected_groups),
            "site_reference_only": len(output_rows) - len(supervised),
            "supervised_s0": len(supervised),
            "unique_supervised_ligands": len({row["inchikey"] for row in supervised}),
        },
        "rejections": dict(sorted(rejections.items())),
        "replicate_audit": replicate_records,
        "output": {
            "path": str(output_tsv.resolve()),
            "sha256": sha256_bytes(payload),
            "bytes": len(payload),
        },
    }
    preserve_manifest_timestamp(manifest_path, manifest, "created_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    return manifest


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def download_coordinates(
    selection_tsv: Path,
    coordinate_root: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    pdb_ids = sorted({row["pdb_id"].upper() for row in read_tsv(selection_tsv)})
    records: list[dict[str, Any]] = []
    coordinate_root.mkdir(parents=True, exist_ok=True)
    for index, pdb_id in enumerate(pdb_ids):
        destination = coordinate_root / f"{pdb_id}.cif.gz"
        url = COORDINATE_URL.format(pdb_id=pdb_id)
        if not destination.exists():
            request = Request(url, headers={"User-Agent": USER_AGENT})
            try:
                with urlopen(request, timeout=120) as response:
                    payload = response.read()
            except (HTTPError, URLError, TimeoutError) as error:
                raise RcsbPilotError(f"Cannot download {pdb_id}: {error}") from error
            try:
                decompressed = gzip.decompress(payload)
            except (OSError, EOFError) as error:
                raise RcsbPilotError(f"Invalid gzip coordinate payload for {pdb_id}") from error
            if not decompressed.startswith(b"data_"):
                raise RcsbPilotError(f"Coordinate payload for {pdb_id} is not mmCIF")
            immutable_write(destination, payload)
            if index and index % 50 == 0:
                time.sleep(0.25)
        records.append(
            {
                "pdb_id": pdb_id,
                "url": url,
                "path": str(destination.resolve()),
                "sha256": sha256_file(destination),
                "bytes": destination.stat().st_size,
            }
        )
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "source": "RCSB PDB coordinate download service",
        "selection_tsv": str(selection_tsv.resolve()),
        "selection_sha256": sha256_file(selection_tsv),
        "coordinate_format": "PDBx/mmCIF gzip",
        "coordinate_count": len(records),
        "records": records,
    }
    preserve_manifest_timestamp(manifest_path, manifest, "created_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    return manifest
