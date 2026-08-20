"""Cross-check the bounded RCSB pilot against an immutable BindingDB TSV archive."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


REQUIRED_BINDINGDB_COLUMNS = (
    "BindingDB Reactant_set_id",
    "Ligand InChI Key",
    "Kd (nM)",
    "Curation/DataSource",
    "Article DOI",
    "PMID",
    "Date of publication",
    "Date in BindingDB",
    "Ligand HET ID in PDB",
    "PDB ID(s) for Ligand-Target Complex",
)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_tsv_line(raw: bytes) -> list[str]:
    decoded = raw.decode("utf-8", errors="replace").rstrip("\r\n")
    return next(csv.reader([decoded], delimiter="\t", quotechar='"'))


def _normalise_doi(value: str) -> str:
    cleaned = value.strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if cleaned.startswith(prefix):
            cleaned = cleaned[len(prefix) :]
    return cleaned.rstrip(".")


def _normalise_pmid(value: str) -> str:
    match = re.search(r"\d+", value or "")
    return match.group(0) if match else ""


def _parse_exact_number(value: str) -> float | None:
    cleaned = value.strip().replace(",", "")
    if not cleaned or cleaned[0] in "<>~":
        return None
    try:
        number = float(cleaned)
    except ValueError:
        return None
    return number if math.isfinite(number) and number > 0 else None


def _pdb_ids(value: str) -> set[str]:
    return {token.upper() for token in re.findall(r"(?i)(?<![a-z0-9])[a-z0-9]{4}(?![a-z0-9])", value)}


def _uniprot_ids(record: dict[str, str]) -> set[str]:
    result: set[str] = set()
    for name, value in record.items():
        if "UniProt" not in name or "Primary ID of Target Chain" not in name:
            continue
        for token in re.split(r"[\s,;|]+", value.strip()):
            if token:
                result.add(token.upper())
    return result


def load_pilot(path: Path) -> list[dict[str, str]]:
    """Load and minimally validate the frozen RCSB pilot."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "pdb_id",
        "ligand_comp_id",
        "value_nm",
        "inchikey",
        "uniprot_accession",
        "citation_doi",
        "citation_pubmed",
        "observation_id",
        "role",
        "label_quarantined",
    }
    missing = required.difference(rows[0] if rows else {})
    if missing:
        raise ValueError(f"pilot is missing columns: {sorted(missing)}")
    if not rows:
        raise ValueError("pilot is empty")
    return rows


def scan_bindingdb(
    archive: Path,
    wanted_pdb_ids: Iterable[str],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Scan one zipped BindingDB TSV and retain only rows mentioning pilot PDB IDs."""
    wanted = {value.upper() for value in wanted_pdb_ids}
    if not wanted:
        raise ValueError("no pilot PDB IDs were supplied")
    retained: list[dict[str, str]] = []
    scanned_rows = 0
    with zipfile.ZipFile(archive) as bundle:
        names = bundle.namelist()
        if len(names) != 1:
            raise ValueError(f"expected one BindingDB TSV in archive, found {names}")
        info = bundle.getinfo(names[0])
        with bundle.open(info) as binary:
            header = _parse_tsv_line(binary.readline())
            missing = set(REQUIRED_BINDINGDB_COLUMNS).difference(header)
            if missing:
                raise ValueError(f"BindingDB archive is missing columns: {sorted(missing)}")
            pdb_column = header.index("PDB ID(s) for Ligand-Target Complex")
            for raw in binary:
                scanned_rows += 1
                prefix = raw.split(b"\t", pdb_column + 1)
                if len(prefix) <= pdb_column or not prefix[pdb_column]:
                    continue
                mentioned = {
                    token.decode("ascii").upper()
                    for token in re.findall(
                        rb"(?i)(?<![a-z0-9])[a-z0-9]{4}(?![a-z0-9])",
                        prefix[pdb_column],
                    )
                }
                if not mentioned.intersection(wanted):
                    continue
                values = _parse_tsv_line(raw)
                if len(values) < len(header):
                    values.extend([""] * (len(header) - len(values)))
                record = dict(zip(header, values))
                if _pdb_ids(record["PDB ID(s) for Ligand-Target Complex"]) & wanted:
                    retained.append(record)

    scan = {
        "archive_entry": names[0],
        "archive_entry_uncompressed_bytes": info.file_size,
        "archive_entry_crc32": f"{info.CRC:08x}",
        "scanned_data_rows": scanned_rows,
        "retained_candidate_rows": len(retained),
    }
    return retained, scan


def _candidate_assessment(pilot: dict[str, str], candidate: dict[str, str]) -> dict[str, Any]:
    expected_kd = _parse_exact_number(pilot["value_nm"])
    observed_kd = _parse_exact_number(candidate["Kd (nM)"])
    pilot_doi = _normalise_doi(pilot["citation_doi"])
    candidate_doi = _normalise_doi(candidate["Article DOI"])
    pilot_pmid = _normalise_pmid(pilot["citation_pubmed"])
    candidate_pmid = _normalise_pmid(candidate["PMID"])

    pdb_match = pilot["pdb_id"].upper() in _pdb_ids(
        candidate["PDB ID(s) for Ligand-Target Complex"]
    )
    het_match = pilot["ligand_comp_id"].upper() in {
        token.upper()
        for token in re.split(r"[\s,;|]+", candidate["Ligand HET ID in PDB"].strip())
        if token
    }
    inchikey_match = (
        bool(pilot["inchikey"])
        and pilot["inchikey"].upper() == candidate["Ligand InChI Key"].strip().upper()
    )
    ligand_match = het_match or inchikey_match
    uniprot_match = pilot["uniprot_accession"].upper() in _uniprot_ids(candidate)
    same_publication = bool(
        (pilot_doi and pilot_doi == candidate_doi)
        or (pilot_pmid and pilot_pmid == candidate_pmid)
    )
    measurement_publication_present = bool(candidate_doi or candidate_pmid)
    kd_match = bool(
        expected_kd is not None
        and observed_kd is not None
        and math.isclose(expected_kd, observed_kd, rel_tol=1e-9, abs_tol=1e-9)
    )
    supervised = pilot["role"] == "supervised_s0"
    high_confidence = bool(
        pdb_match
        and ligand_match
        and uniprot_match
        and measurement_publication_present
        and (kd_match if supervised else True)
    )
    score = (
        3 * sum((pdb_match, ligand_match, uniprot_match, kd_match))
        + int(measurement_publication_present)
        + int(same_publication)
    )
    return {
        "candidate": candidate,
        "score": score,
        "pdb_match": pdb_match,
        "het_match": het_match,
        "inchikey_match": inchikey_match,
        "ligand_match": ligand_match,
        "uniprot_match": uniprot_match,
        "measurement_publication_present": measurement_publication_present,
        "same_structure_measurement_publication": same_publication,
        "kd_match": kd_match,
        "high_confidence": high_confidence,
    }


def reconcile(
    pilot_rows: list[dict[str, str]],
    candidates: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Select the best BindingDB provenance row for each pilot observation."""
    by_pdb: dict[str, list[dict[str, str]]] = {}
    for candidate in candidates:
        for pdb_id in _pdb_ids(candidate["PDB ID(s) for Ligand-Target Complex"]):
            by_pdb.setdefault(pdb_id, []).append(candidate)

    output: list[dict[str, Any]] = []
    source_counts: Counter[str] = Counter()
    for pilot in pilot_rows:
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
        best = assessments[0] if assessments else None
        tied = sum(item["score"] == best["score"] for item in assessments) if best else 0
        candidate = best["candidate"] if best else {}
        if candidate:
            source_counts[candidate["Curation/DataSource"] or "<blank>"] += 1
        supervised = pilot["role"] == "supervised_s0"
        output.append(
            {
                "observation_id": pilot["observation_id"],
                "role": pilot["role"],
                "pdb_id": pilot["pdb_id"],
                "ligand_comp_id": pilot["ligand_comp_id"],
                "bindingdb_reactant_set_id": candidate.get("BindingDB Reactant_set_id", ""),
                "candidate_count_for_pdb": len(assessments),
                "top_score_tie_count": tied,
                "pdb_match": bool(best and best["pdb_match"]),
                "ligand_het_match": bool(best and best["het_match"]),
                "ligand_inchikey_match": bool(best and best["inchikey_match"]),
                "uniprot_match": bool(best and best["uniprot_match"]),
                "measurement_publication_present": bool(
                    best and best["measurement_publication_present"]
                ),
                "same_structure_measurement_publication": bool(
                    best and best["same_structure_measurement_publication"]
                ),
                "kd_exact_match": bool(best and best["kd_match"]) if supervised else "",
                "bindingdb_kd_status": (
                    "exact"
                    if supervised and best and best["kd_match"]
                    else (
                        "censored_or_missing"
                        if supervised
                        else "quarantined_reference_not_exported"
                    )
                ),
                "high_confidence_provenance": bool(best and best["high_confidence"]),
                "curation_data_source": candidate.get("Curation/DataSource", ""),
                "article_doi": candidate.get("Article DOI", ""),
                "pmid": candidate.get("PMID", ""),
                "publication_date": candidate.get("Date of publication", ""),
                "bindingdb_curation_date": candidate.get("Date in BindingDB", ""),
                "label_quarantine_preserved": (
                    pilot["label_quarantined"].strip().lower() == "true"
                    and pilot["value_nm"] == ""
                    and pilot["role"] != "supervised_s0"
                )
                if pilot["role"] != "supervised_s0"
                else "",
            }
        )

    supervised_rows = [row for row in output if row["role"] == "supervised_s0"]
    reference_rows = [row for row in output if row["role"] != "supervised_s0"]
    summary = {
        "pilot_rows": len(output),
        "supervised_rows": len(supervised_rows),
        "reference_only_rows": len(reference_rows),
        "supervised_high_confidence": sum(
            bool(row["high_confidence_provenance"]) for row in supervised_rows
        ),
        "supervised_exact_kd": sum(row["kd_exact_match"] is True for row in supervised_rows),
        "supervised_ligand_identity": sum(
            bool(row["ligand_het_match"] or row["ligand_inchikey_match"])
            for row in supervised_rows
        ),
        "supervised_uniprot": sum(bool(row["uniprot_match"]) for row in supervised_rows),
        "supervised_measurement_publication_present": sum(
            bool(row["measurement_publication_present"]) for row in supervised_rows
        ),
        "same_structure_measurement_publication": sum(
            bool(row["same_structure_measurement_publication"]) for row in supervised_rows
        ),
        "reference_label_quarantine_preserved": sum(
            row["label_quarantine_preserved"] is True for row in reference_rows
        ),
        "ambiguous_top_score_rows": sum(row["top_score_tie_count"] > 1 for row in output),
        "curation_data_sources": dict(sorted(source_counts.items())),
    }
    return output, summary


def build_strict_dataset(
    pilot_rows: list[dict[str, str]],
    provenance_rows: list[dict[str, Any]],
    pocket_path: Path,
    strict_pilot_path: Path,
    strict_pocket_path: Path,
    min_supervised_per_construct: int,
) -> dict[str, Any]:
    """Create a label-safe view after provenance and minimum-group filtering."""
    provenance_by_id = {str(row["observation_id"]): row for row in provenance_rows}
    eligible = {
        observation_id
        for observation_id, row in provenance_by_id.items()
        if row["role"] == "supervised_s0" and row["high_confidence_provenance"] is True
    }
    group_counts = Counter(
        row["construct_group_id"]
        for row in pilot_rows
        if row["observation_id"] in eligible
    )
    retained_groups = {
        group for group, count in group_counts.items() if count >= min_supervised_per_construct
    }
    retained_ids = {
        row["observation_id"]
        for row in pilot_rows
        if row["observation_id"] in eligible and row["construct_group_id"] in retained_groups
    }

    strict_rows: list[dict[str, Any]] = []
    for row in pilot_rows:
        is_supervised = row["role"] == "supervised_s0"
        if is_supervised and row["observation_id"] not in retained_ids:
            continue
        if not is_supervised and row["construct_group_id"] not in retained_groups:
            continue
        output = dict(row)
        if is_supervised:
            provenance = provenance_by_id[row["observation_id"]]
            output.update(
                {
                    "bindingdb_release": "202608",
                    "bindingdb_reactant_set_id": provenance["bindingdb_reactant_set_id"],
                    "measurement_source": provenance["curation_data_source"],
                    "measurement_publication_doi": provenance["article_doi"],
                    "measurement_publication_pmid": provenance["pmid"],
                    "measurement_publication_date": provenance["publication_date"],
                    "bindingdb_curation_date": provenance["bindingdb_curation_date"],
                    "measurement_provenance_status": "exact_pdb_ligand_target_kd_and_publication",
                }
            )
        else:
            output.update(
                {
                    "bindingdb_release": "",
                    "bindingdb_reactant_set_id": "",
                    "measurement_source": "",
                    "measurement_publication_doi": "",
                    "measurement_publication_pmid": "",
                    "measurement_publication_date": "",
                    "bindingdb_curation_date": "",
                    "measurement_provenance_status": "label_quarantined_site_reference_only",
                }
            )
        strict_rows.append(output)
    write_tsv(strict_rows, strict_pilot_path)

    pocket_rows = [
        row for row in load_table(pocket_path) if row["observation_id"] in retained_ids
    ]
    write_tsv(pocket_rows, strict_pocket_path)
    supervised = [row for row in strict_rows if row["role"] == "supervised_s0"]
    references = [row for row in strict_rows if row["role"] != "supervised_s0"]
    if any(row["value_nm"] or row["pKd"] for row in references):
        raise ValueError("strict dataset exposed a site-reference label")
    if len(pocket_rows) != 2 * len(supervised):
        raise ValueError("strict pocket table does not contain exactly S0 and S1 per observation")
    return {
        "status": "PASS",
        "minimum_supervised_per_construct": min_supervised_per_construct,
        "pre_group_filter_high_confidence": len(eligible),
        "retained_supervised": len(supervised),
        "retained_reference_only": len(references),
        "retained_construct_groups": len(retained_groups),
        "retained_group_counts": dict(
            sorted(Counter(row["construct_group_id"] for row in supervised).items())
        ),
        "dropped_low_count_groups": dict(
            sorted(
                (group, count)
                for group, count in group_counts.items()
                if group not in retained_groups
            )
        ),
        "removed_provenance_failures": sum(
            row["role"] == "supervised_s0" and row["observation_id"] not in eligible
            for row in pilot_rows
        ),
        "pilot_output": {
            "path": str(strict_pilot_path),
            "bytes": strict_pilot_path.stat().st_size,
            "sha256": sha256_file(strict_pilot_path),
        },
        "pocket_output": {
            "path": str(strict_pocket_path),
            "bytes": strict_pocket_path.stat().st_size,
            "sha256": sha256_file(strict_pocket_path),
        },
    }


def load_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_audit(
    pilot_path: Path,
    archive_path: Path,
    output_path: Path,
    report_path: Path,
    source_url: str,
    expected_archive_sha256: str = "",
    pocket_path: Path | None = None,
    strict_pilot_path: Path | None = None,
    strict_pocket_path: Path | None = None,
    min_supervised_per_construct: int = 8,
) -> dict[str, Any]:
    """Run, persist, and summarize the complete provenance audit."""
    archive_sha256 = sha256_file(archive_path)
    if expected_archive_sha256 and archive_sha256 != expected_archive_sha256.lower():
        raise ValueError(
            f"BindingDB archive hash mismatch: {archive_sha256} != {expected_archive_sha256}"
        )
    pilot_rows = load_pilot(pilot_path)
    candidates, scan = scan_bindingdb(archive_path, (row["pdb_id"] for row in pilot_rows))
    rows, reconciliation = reconcile(pilot_rows, candidates)
    write_tsv(rows, output_path)
    strict_dataset = None
    strict_args = (pocket_path, strict_pilot_path, strict_pocket_path)
    if any(value is not None for value in strict_args):
        if any(value is None for value in strict_args):
            raise ValueError("all strict-dataset paths must be supplied together")
        strict_dataset = build_strict_dataset(
            pilot_rows,
            rows,
            pocket_path=pocket_path,
            strict_pilot_path=strict_pilot_path,
            strict_pocket_path=strict_pocket_path,
            min_supervised_per_construct=min_supervised_per_construct,
        )

    report = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "PASS"
            if reconciliation["supervised_high_confidence"]
            == reconciliation["supervised_rows"]
            and reconciliation["reference_label_quarantine_preserved"]
            == reconciliation["reference_only_rows"]
            else "FAIL"
        ),
        "source": {
            "name": "BindingDB",
            "release": "202608",
            "url": source_url,
            "archive_path": str(archive_path),
            "archive_bytes": archive_path.stat().st_size,
            "archive_sha256": archive_sha256,
        },
        "pilot": {
            "path": str(pilot_path),
            "sha256": sha256_file(pilot_path),
        },
        "scan": scan,
        "reconciliation": reconciliation,
        "strict_dataset": strict_dataset,
        "output": {
            "path": str(output_path),
            "bytes": output_path.stat().st_size,
            "sha256": sha256_file(output_path),
        },
        "interpretation": (
            "A pass proves record-level agreement between the RCSB pilot annotation and the "
            "pinned BindingDB release on PDB identifier, ligand identity, target accession, "
            "publication, and exact Kd. It does not independently validate the experimental "
            "measurement reported by the source article. Reference-only labels remain quarantined."
        ),
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report
