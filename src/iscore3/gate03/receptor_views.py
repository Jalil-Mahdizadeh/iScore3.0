"""Coordinate-safe S1/S2 receptor views for the frozen Gate-3 cohort.

The historical ligand is represented only by the already-frozen site residue
indices.  This module never accepts a query ligand coordinate or affinity label.
It keeps target-sequence, PDB-entity, and predicted-model residue coordinates
explicitly separate.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np

from iscore3.data.rcsb_gate01 import (
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
    utc_now,
)
from iscore3.protein.pocket_features import (
    AA3_TO_1,
    FEATURE_NAMES,
    AtomRecord,
    PocketInstance,
    extract_protein_pocket,
    pocket_feature_dict,
    read_mmcif_atoms,
    tsv_bytes,
)
from iscore3.protein.structure_views import (
    AA1_TO_3,
    align_construct_to_prediction,
    write_ca_pdb,
)


USER_AGENT = "iScore3.0-gate03/1.0 (scientific provenance audit)"
AF_ENDPOINT = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

CHEMISTRY_GROUPS = {
    "hydrophobic": set("AVILMFWY"),
    "aromatic": set("FWYH"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "polar_uncharged": set("STNQCY"),
    "small": set("AGSTC"),
    "gly_pro": set("GP"),
}
RADIAL_SHELLS = ((0.0, 0.35), (0.35, 0.65), (0.65, float("inf")))
CONTACT_THRESHOLDS = (6.0, 8.0, 10.0, 12.0)

POCKET_V2_NAMES = (
    *FEATURE_NAMES,
    *(f"chem_fraction_{name}" for name in CHEMISTRY_GROUPS),
    *(
        f"radial_shell_{shell}_{name}_fraction"
        for shell in range(len(RADIAL_SHELLS))
        for name in ("all", *CHEMISTRY_GROUPS)
    ),
    *(f"centroid_distance_quantile_{value}" for value in (10, 25, 50, 75, 90)),
    *(
        f"contact_{stat}_{int(threshold)}A"
        for threshold in CONTACT_THRESHOLDS
        for stat in ("density", "mean_degree", "degree_std", "maximum_degree")
    ),
    "shape_linearity",
    "shape_planarity",
    "shape_sphericity",
    "shape_anisotropy",
    "shape_convex_hull_volume",
    "shape_residues_per_hull_volume",
    "exposure_low_contact_fraction_8A",
    "exposure_high_contact_fraction_8A",
    "quality_ca_bfactor_mean",
    "quality_ca_bfactor_std",
    "quality_ca_bfactor_min",
    "quality_ca_bfactor_q10",
)


class Gate3ViewError(RuntimeError):
    """Raised when a Gate-3 receptor-view invariant fails."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _download(path: Path, url: str) -> bytes:
    if path.exists():
        return path.read_bytes()
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=180) as response:
        payload = response.read()
    immutable_write(path, payload)
    return payload


def _accessions(value: str) -> tuple[str, ...]:
    accessions = sorted({part.strip() for part in value.split(";") if part.strip()})
    if not accessions:
        raise Gate3ViewError("No UniProt accession declared")
    return tuple(accessions)


def acquire_alphafold(
    dataset: Path,
    raw_root: Path,
    manifest_path: Path,
    endpoint_template: str = AF_ENDPOINT,
) -> dict[str, Any]:
    """Acquire one exact canonical AlphaFold DB model per Gate-3 series."""

    rows = read_tsv(dataset)
    accessions = sorted(
        {accession for row in rows for accession in _accessions(row["uniprot_ids"])}
    )
    records: list[dict[str, Any]] = []
    for accession in accessions:
        api_url = endpoint_template.format(accession=accession)
        api_path = raw_root / "api" / f"{accession}.json"
        try:
            api_payload = _download(api_path, api_url)
            response = json.loads(api_payload)
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as error:
            records.append(
                {
                    "uniprot_accession": accession,
                    "status": "api_unavailable",
                    "api_url": api_url,
                    "error": f"{type(error).__name__}: {error}",
                }
            )
            continue
        exact = [
            item
            for item in response
            if item.get("uniprotAccession") == accession
            and item.get("isUniProt") is True
            and item.get("isComplex") is not True
        ]
        if not exact:
            records.append(
                {
                    "uniprot_accession": accession,
                    "status": "exact_canonical_accession_unavailable",
                    "api_url": api_url,
                    "api_path": str(api_path),
                    "api_sha256": hashlib.sha256(api_payload).hexdigest(),
                }
            )
            continue
        selected = max(
            exact,
            key=lambda item: (
                int(item.get("latestVersion") or -1),
                str(item.get("modelCreatedDate") or ""),
                str(item.get("entryId") or ""),
            ),
        )
        sequence = str(selected.get("sequence") or selected.get("uniprotSequence") or "")
        cif_url = str(selected.get("cifUrl") or "")
        if not sequence or not cif_url:
            raise Gate3ViewError(f"Incomplete AlphaFold record for {accession}")
        cif_path = raw_root / "structures" / Path(cif_url).name
        cif_payload = _download(cif_path, cif_url)
        records.append(
            {
                "uniprot_accession": accession,
                "status": "available",
                "api_url": api_url,
                "api_path": str(api_path),
                "api_bytes": len(api_payload),
                "api_sha256": hashlib.sha256(api_payload).hexdigest(),
                "entry_id": selected.get("entryId"),
                "latest_version": selected.get("latestVersion"),
                "model_created_date": selected.get("modelCreatedDate"),
                "tool_used": selected.get("toolUsed"),
                "global_metric_value": selected.get("globalMetricValue"),
                "sequence": sequence,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "sequence_length": len(sequence),
                "cif_url": cif_url,
                "cif_path": str(cif_path),
                "cif_bytes": len(cif_payload),
                "cif_sha256": hashlib.sha256(cif_payload).hexdigest(),
            }
        )
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "source": "AlphaFold Protein Structure Database API",
        "endpoint_template": endpoint_template,
        "snapshot_semantics": "mutable_API_response_and_versioned_model_cached_immutably",
        "dataset": {"path": str(dataset), "sha256": sha256_file(dataset)},
        "counts": {
            "requested_accessions": len(accessions),
            "available": sum(row["status"] == "available" for row in records),
            "unavailable": sum(row["status"] != "available" for row in records),
        },
        "records": records,
    }
    preserve_manifest_timestamp(manifest_path, manifest, "created_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    return manifest


def _residue_centroids(pocket: PocketInstance) -> tuple[list[str], np.ndarray]:
    names = []
    centroids = []
    for position in pocket.present_positions:
        residue_atoms = [atom for atom in pocket.atoms if atom.seq_id == position]
        names.append(AA3_TO_1.get(pocket.residue_name_by_position[position], "X"))
        centroids.append(np.mean([atom.xyz for atom in residue_atoms], axis=0))
    return names, np.asarray(centroids, dtype=np.float64)


def pocket_descriptor_v2(pocket: PocketInstance) -> dict[str, float]:
    """Rigid-invariant pocket geometry, chemistry, topology, and quality features."""

    values = dict(pocket_feature_dict(pocket))
    names, centroids = _residue_centroids(pocket)
    count = len(names)
    chemistry = {
        name: np.asarray([residue in residues for residue in names], dtype=bool)
        for name, residues in CHEMISTRY_GROUPS.items()
    }
    for name, mask in chemistry.items():
        values[f"chem_fraction_{name}"] = float(np.mean(mask))

    centre = np.mean(centroids, axis=0)
    radii = np.linalg.norm(centroids - centre, axis=1)
    scale = max(float(np.max(radii)), 1.0e-12)
    normalized = radii / scale
    for shell_index, (lower, upper) in enumerate(RADIAL_SHELLS):
        shell = (normalized >= lower) & (normalized < upper)
        values[f"radial_shell_{shell_index}_all_fraction"] = float(np.mean(shell))
        for name, mask in chemistry.items():
            values[f"radial_shell_{shell_index}_{name}_fraction"] = float(
                np.mean(shell & mask)
            )

    distances = np.sqrt(
        np.sum((centroids[:, None, :] - centroids[None, :, :]) ** 2, axis=2)
    )
    distances = np.round(distances, decimals=8)
    upper = distances[np.triu_indices(count, k=1)]
    for quantile in (10, 25, 50, 75, 90):
        values[f"centroid_distance_quantile_{quantile}"] = float(
            np.percentile(upper, quantile)
        )
    degrees_8a = None
    for threshold in CONTACT_THRESHOLDS:
        adjacency = (distances <= threshold) & (~np.eye(count, dtype=bool))
        degrees = adjacency.sum(axis=1).astype(np.float64)
        prefix = f"contact_{{}}_{int(threshold)}A"
        values[prefix.format("density")] = float(adjacency.sum() / (count * (count - 1)))
        values[prefix.format("mean_degree")] = float(np.mean(degrees))
        values[prefix.format("degree_std")] = float(np.std(degrees))
        values[prefix.format("maximum_degree")] = float(np.max(degrees))
        if threshold == 8.0:
            degrees_8a = degrees
    assert degrees_8a is not None

    covariance = np.cov(centroids, rowvar=False, bias=True)
    eigenvalues = np.clip(np.linalg.eigvalsh(covariance)[::-1], 0.0, None)
    largest = max(float(eigenvalues[0]), 1.0e-12)
    values["shape_linearity"] = float((eigenvalues[0] - eigenvalues[1]) / largest)
    values["shape_planarity"] = float((eigenvalues[1] - eigenvalues[2]) / largest)
    values["shape_sphericity"] = float(eigenvalues[2] / largest)
    values["shape_anisotropy"] = float(1.0 - eigenvalues[2] / largest)
    try:
        from scipy.spatial import ConvexHull

        hull_volume = float(ConvexHull(centroids).volume) if count >= 4 else 0.0
    except Exception:
        hull_volume = 0.0
    values["shape_convex_hull_volume"] = hull_volume
    values["shape_residues_per_hull_volume"] = float(count / max(hull_volume, 1.0))
    values["exposure_low_contact_fraction_8A"] = float(np.mean(degrees_8a <= 2.0))
    values["exposure_high_contact_fraction_8A"] = float(
        np.mean(degrees_8a >= np.median(degrees_8a))
    )

    ca_b = np.asarray(
        [
            atom.b_factor
            for atom in pocket.atoms
            if atom.atom_name.strip() == "CA" and math.isfinite(atom.b_factor)
        ],
        dtype=np.float64,
    )
    if len(ca_b):
        quality = (np.mean(ca_b), np.std(ca_b), np.min(ca_b), np.percentile(ca_b, 10))
    else:
        quality = (0.0, 0.0, 0.0, 0.0)
    for name, value in zip(
        (
            "quality_ca_bfactor_mean",
            "quality_ca_bfactor_std",
            "quality_ca_bfactor_min",
            "quality_ca_bfactor_q10",
        ),
        quality,
        strict=True,
    ):
        values[name] = float(value)

    if tuple(values) != POCKET_V2_NAMES or not all(math.isfinite(v) for v in values.values()):
        raise Gate3ViewError("Pocket-v2 feature schema/order or finiteness failure")
    return values


def _protein_chain(atoms: Sequence[AtomRecord]) -> tuple[str, str, list[AtomRecord]]:
    candidates: dict[tuple[str, str], list[AtomRecord]] = defaultdict(list)
    for atom in atoms:
        if atom.group == "ATOM" and atom.seq_id is not None:
            candidates[(atom.entity_id, atom.asym_id)].append(atom)
    if not candidates:
        raise Gate3ViewError("No protein chain in predicted structure")
    entity, asym = min(
        candidates,
        key=lambda key: (
            -len({atom.seq_id for atom in candidates[key] if atom.atom_name.strip() == "CA"}),
            key,
        ),
    )
    return entity, asym, candidates[(entity, asym)]


def _site_definition(
    definitions: Sequence[Mapping[str, Any]], series_id: str, pdb_id: str
) -> Mapping[str, Any]:
    matches = [
        row
        for row in definitions
        if row["construct_group_id"] == series_id and row["reference_pdb_id"] == pdb_id
    ]
    if len(matches) != 1:
        raise Gate3ViewError(
            f"Expected one frozen site for {series_id}/{pdb_id}, found {len(matches)}"
        )
    return matches[0]


def _feature_row(
    *,
    series_id: str,
    view: str,
    source_id: str,
    source_sha256: str,
    entity_id: str,
    asym_id: str,
    pocket: PocketInstance,
    target_positions: Sequence[int],
    mean_site_plddt: float | None,
) -> dict[str, Any]:
    return {
        "series_id": series_id,
        "view": view,
        "source_structure_id": source_id,
        "source_structure_sha256": source_sha256,
        "protein_entity_id": entity_id,
        "protein_asym_id": asym_id,
        "source_site_positions": ";".join(map(str, pocket.expected_positions)),
        "target_site_positions": ";".join(map(str, target_positions)),
        "present_source_positions": ";".join(map(str, pocket.present_positions)),
        "site_mean_plddt": "" if mean_site_plddt is None else mean_site_plddt,
        "query_ligand_coordinates_read": False,
        "affinity_labels_read": False,
        **pocket_descriptor_v2(pocket),
    }


def build_s1_s2_views(
    *,
    dataset: Path,
    sites: Path,
    experimental_coordinate_root: Path,
    alphafold_manifest: Path,
    derived_root: Path,
    feature_output: Path,
    view_manifest: Path,
    audit_output: Path,
    minimum_site_coverage: float = 0.80,
    minimum_site_plddt: float = 70.0,
) -> dict[str, Any]:
    """Build one S1 and eligible S2 view per frozen assay series."""

    rows = read_tsv(dataset)
    by_series: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_series[row["series_id"]].append(row)
    definitions = json.loads(sites.read_text(encoding="utf-8"))["definitions"]
    af = {
        row["uniprot_accession"]: row
        for row in json.loads(alphafold_manifest.read_text(encoding="utf-8"))["records"]
        if row["status"] == "available"
    }
    feature_rows: list[dict[str, Any]] = []
    views: list[dict[str, Any]] = []
    for series_id in sorted(by_series):
        representative = by_series[series_id][0]
        sequence = representative["target_sequence"]
        reference_id = representative["site_reference_pdb_id"]
        site = _site_definition(definitions, series_id, reference_id)
        reference_path = experimental_coordinate_root / f"{reference_id}.cif.gz"
        if sha256_file(reference_path) != representative["site_reference_structure_sha256"]:
            raise Gate3ViewError(f"Reference coordinate hash mismatch for {series_id}")
        reference_atoms = read_mmcif_atoms(reference_path)
        expected_names = {
            int(key): value for key, value in site["residue_name_by_position"].items()
        }
        s1 = extract_protein_pocket(
            reference_atoms,
            pdb_id=reference_id,
            structure_path=reference_path,
            protein_entity_id=site["reference_protein_entity_id"],
            candidate_asym_ids=[site["reference_protein_asym_id"]],
            positions=site["positions_label_seq_id"],
            expected_residue_names=expected_names,
            minimum_coverage=minimum_site_coverage,
        )
        target_positions = [int(value) for value in site["target_positions"]]
        if target_positions != [int(value) for value in representative["site_target_positions"].split(";")]:
            raise Gate3ViewError(f"Target-site ledger mismatch for {series_id}")
        feature_rows.append(
            _feature_row(
                series_id=series_id,
                view="S1",
                source_id=reference_id,
                source_sha256=sha256_file(reference_path),
                entity_id=site["reference_protein_entity_id"],
                asym_id=s1.selected_asym_id,
                pocket=s1,
                target_positions=target_positions,
                mean_site_plddt=None,
            )
        )
        s1_chain = [
            atom
            for atom in reference_atoms
            if atom.group == "ATOM"
            and atom.entity_id == site["reference_protein_entity_id"]
            and atom.asym_id == s1.selected_asym_id
        ]
        s1_global = write_ca_pdb(s1_chain, derived_root / "S1" / f"{series_id}-global.pdb")
        s1_local = write_ca_pdb(s1.atoms, derived_root / "S1" / f"{series_id}-pocket.pdb")
        views.append(
            {
                "series_id": series_id,
                "view": "S1",
                "status": "eligible",
                "source_id": reference_id,
                "source_structure_path": str(reference_path),
                "source_structure_sha256": sha256_file(reference_path),
                "protein_entity_id": site["reference_protein_entity_id"],
                "protein_asym_id": s1.selected_asym_id,
                "source_site_positions": list(s1.expected_positions),
                "target_site_positions": target_positions,
                "global": s1_global,
                "pocket": s1_local,
            }
        )

        candidates = [
            af[accession]
            for accession in _accessions(representative["uniprot_ids"])
            if accession in af
        ]
        scored_candidates = []
        for candidate in candidates:
            candidate_alignment = align_construct_to_prediction(
                sequence, candidate["sequence"]
            )
            scored_candidates.append(
                (candidate_alignment["aligned_identity"],
                 candidate_alignment["construct_coverage"],
                 candidate["uniprot_accession"], candidate, candidate_alignment)
            )
        selected = max(scored_candidates, default=None)
        prediction = selected[3] if selected is not None else None
        if prediction is None:
            views.append({"series_id": series_id, "view": "S2", "status": "unavailable"})
            continue
        alignment = selected[4]
        site_map = {
            position: alignment["mapping"][position]
            for position in target_positions
            if position in alignment["mapping"]
        }
        site_coverage = len(site_map) / len(target_positions)
        site_identity = (
            sum(
                sequence[source - 1] == prediction["sequence"][target - 1]
                for source, target in site_map.items()
            )
            / len(site_map)
            if site_map
            else 0.0
        )
        mapping = {
            key: value for key, value in alignment.items() if key != "mapping"
        }
        mapping.update(
            {
                "declared_uniprot_accessions": list(
                    _accessions(representative["uniprot_ids"])
                ),
                "selected_uniprot_accession": prediction["uniprot_accession"],
                "site_mapping_coverage": site_coverage,
                "site_identity": site_identity,
                "target_to_prediction_site_positions": {
                    str(key): value for key, value in sorted(site_map.items())
                },
            }
        )
        if (
            alignment["construct_coverage"] < 0.80
            or alignment["aligned_identity"] < 0.98
            or site_coverage < minimum_site_coverage
            or site_identity < 0.98
        ):
            views.append(
                {
                    "series_id": series_id,
                    "view": "S2",
                    "status": "sequence_mapping_threshold_failure",
                    "mapping": mapping,
                }
            )
            continue
        af_path = Path(prediction["cif_path"])
        if sha256_file(af_path) != prediction["cif_sha256"]:
            raise Gate3ViewError(f"AlphaFold coordinate hash mismatch for {series_id}")
        af_atoms = read_mmcif_atoms(af_path)
        entity_id, asym_id, chain_atoms = _protein_chain(af_atoms)
        predicted_positions = sorted(site_map.values())
        predicted_names = {
            position: AA1_TO_3.get(prediction["sequence"][position - 1], "UNK")
            for position in predicted_positions
        }
        try:
            s2 = extract_protein_pocket(
                af_atoms,
                pdb_id=str(prediction["entry_id"]),
                structure_path=af_path,
                protein_entity_id=entity_id,
                candidate_asym_ids=[asym_id],
                positions=predicted_positions,
                expected_residue_names=predicted_names,
                minimum_coverage=minimum_site_coverage,
            )
        except Exception as error:
            views.append(
                {
                    "series_id": series_id,
                    "view": "S2",
                    "status": "coordinate_mapping_failure",
                    "mapping": mapping,
                    "error": str(error),
                }
            )
            continue
        site_ca = [
            atom.b_factor
            for atom in s2.atoms
            if atom.atom_name.strip() == "CA" and math.isfinite(atom.b_factor)
        ]
        mean_plddt = float(np.mean(site_ca)) if site_ca else float("nan")
        if not math.isfinite(mean_plddt) or mean_plddt < minimum_site_plddt:
            views.append(
                {
                    "series_id": series_id,
                    "view": "S2",
                    "status": "site_plddt_threshold_failure",
                    "mapping": mapping,
                    "mean_site_plddt": mean_plddt,
                }
            )
            continue
        feature_rows.append(
            _feature_row(
                series_id=series_id,
                view="S2",
                source_id=str(prediction["entry_id"]),
                source_sha256=sha256_file(af_path),
                entity_id=entity_id,
                asym_id=asym_id,
                pocket=s2,
                target_positions=target_positions,
                mean_site_plddt=mean_plddt,
            )
        )
        s2_global = write_ca_pdb(chain_atoms, derived_root / "S2" / f"{series_id}-global.pdb")
        s2_local = write_ca_pdb(s2.atoms, derived_root / "S2" / f"{series_id}-pocket.pdb")
        views.append(
            {
                "series_id": series_id,
                "view": "S2",
                "status": "eligible",
                "source_id": prediction["entry_id"],
                "source_structure_path": str(af_path),
                "source_structure_sha256": sha256_file(af_path),
                "protein_entity_id": entity_id,
                "protein_asym_id": asym_id,
                "source_site_positions": list(s2.expected_positions),
                "target_site_positions": target_positions,
                "mean_site_plddt": mean_plddt,
                "mapping": mapping,
                "global": s2_global,
                "pocket": s2_local,
            }
        )

    feature_rows.sort(key=lambda row: (row["series_id"], row["view"]))
    immutable_write(feature_output, tsv_bytes(feature_rows))
    status_counts = Counter((row["view"], row["status"]) for row in views)
    eligible_s2 = status_counts[("S2", "eligible")]
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "information_boundary": {
            "query_ligand_coordinates_read": False,
            "affinity_labels_read": False,
            "site_source": "frozen historical-reference residue indices only",
        },
        "inputs": {
            "dataset": {"path": str(dataset), "sha256": sha256_file(dataset)},
            "sites": {"path": str(sites), "sha256": sha256_file(sites)},
            "alphafold": {
                "path": str(alphafold_manifest),
                "sha256": sha256_file(alphafold_manifest),
            },
        },
        "thresholds": {
            "minimum_target_alignment_coverage": 0.80,
            "minimum_target_aligned_identity": 0.98,
            "minimum_site_mapping_coverage": minimum_site_coverage,
            "minimum_site_identity": 0.98,
            "minimum_mean_site_plddt": minimum_site_plddt,
        },
        "feature_schema": list(POCKET_V2_NAMES),
        "feature_output": {
            "path": str(feature_output),
            "sha256": sha256_file(feature_output),
            "bytes": feature_output.stat().st_size,
        },
        "counts": {
            "series": len(by_series),
            "S1_eligible": status_counts[("S1", "eligible")],
            "S2_eligible": eligible_s2,
            "S2_coverage_fraction": eligible_s2 / len(by_series),
        },
        "views": views,
    }
    preserve_manifest_timestamp(view_manifest, manifest, "created_utc")
    immutable_write(view_manifest, stable_json_bytes(manifest))
    checks = {
        "all_series_have_S1": status_counts[("S1", "eligible")] == len(by_series),
        "S2_series_coverage_at_least_0p80": eligible_s2 / len(by_series) >= 0.80,
        "all_feature_rows_ligand_coordinate_free": all(
            row["query_ligand_coordinates_read"] is False for row in feature_rows
        ),
        "all_feature_rows_label_free": all(row["affinity_labels_read"] is False for row in feature_rows),
        "feature_schema_exact": all(
            all(name in row for name in POCKET_V2_NAMES) for row in feature_rows
        ),
    }
    audit = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": manifest["counts"],
        "view_manifest_sha256": sha256_file(view_manifest),
        "feature_output_sha256": sha256_file(feature_output),
    }
    preserve_manifest_timestamp(audit_output, audit, "created_utc")
    immutable_write(audit_output, stable_json_bytes(audit))
    return audit
