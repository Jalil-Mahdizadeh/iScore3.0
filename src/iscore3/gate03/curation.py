"""BindingDB-first, provenance-preserving Gate-3 assay-series curation."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
from io import StringIO
import math
from pathlib import Path
import re
from statistics import median
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from iscore3.data.rcsb_gate01 import sha256_file


EXACT_NUMBER = re.compile(
    r"^\s*(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\s*$"
)
NUMBER_IN_CONDITION = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)")

REQUIRED_COLUMNS = {
    "BindingDB Reactant_set_id",
    "Ligand SMILES",
    "Ligand InChI Key",
    "Target Name",
    "Kd (nM)",
    "pH",
    "Temp (C)",
    "Curation/DataSource",
    "Article DOI",
    "PMID",
    "Authors",
    "Date of publication",
    "PDB ID(s) for Ligand-Target Complex",
    "Number of Protein Chains in Target (>1 implies a multichain complex)",
    "BindingDB Target Chain Sequence 1",
    "PDB ID(s) of Target Chain 1",
    "UniProt (SwissProt) Primary ID of Target Chain 1",
    "UniProt (TrEMBL) Primary ID of Target Chain 1",
}


class CurationError(RuntimeError):
    """Raised when a Gate-3 curation contract is violated."""


def normalize_doi(value: str) -> str:
    cleaned = (value or "").strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if cleaned.startswith(prefix):
            cleaned = cleaned[len(prefix) :]
    return cleaned.rstrip(" .")


def publication_id(doi: str, pmid: str) -> str:
    normalized_doi = normalize_doi(doi)
    if normalized_doi:
        return f"doi:{normalized_doi}"
    normalized_pmid = (pmid or "").strip()
    if normalized_pmid and normalized_pmid.isdigit():
        return f"pmid:{normalized_pmid}"
    return ""


def normalize_condition(value: str) -> str:
    cleaned = " ".join((value or "").strip().lower().split())
    if not cleaned:
        return "<missing>"
    numbers = NUMBER_IN_CONDITION.findall(cleaned)
    if len(numbers) == 1:
        return f"{float(numbers[0]):.4g}"
    return cleaned


def exact_positive_number(value: str) -> float | None:
    if not EXACT_NUMBER.fullmatch(value or ""):
        return None
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        return None
    return number


def pkd_from_nm(value_nm: float) -> float:
    return 9.0 - math.log10(value_nm)


def split_pdb_ids(value: str) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                token.upper()
                for token in re.findall(r"\b[0-9][A-Za-z0-9]{3}\b", value or "")
            }
        )
    )


def split_accessions(*values: str) -> tuple[str, ...]:
    result = set()
    for value in values:
        result.update(
            token.strip()
            for token in re.split(r"[,;|\s]+", value or "")
            if token.strip()
        )
    return tuple(sorted(result))


def series_identifier(
    sequence_sha256: str, publication: str, ph: str, temperature: str
) -> str:
    payload = "\0".join((sequence_sha256, publication, ph, temperature))
    return "G3S-" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def _canonical_ligand(smiles: str) -> tuple[str, str]:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or molecule.GetNumAtoms() == 0 or molecule.GetNumConformers() != 0:
        raise ValueError("SMILES is not a coordinate-free valid molecule")
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    scaffold_smiles = Chem.MolToSmiles(
        scaffold, canonical=True, isomericSmiles=False
    )
    return canonical, scaffold_smiles


def serialize_tsv(rows: Sequence[Mapping[str, Any]]) -> bytes:
    if not rows:
        raise CurationError("Cannot serialize an empty Gate-3 table")
    fieldnames = list(rows[0])
    existing = set(fieldnames)
    fieldnames.extend(
        sorted({key for row in rows for key in row}.difference(existing))
    )
    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def _percentile(values: Sequence[float], quantile: float) -> float:
    return float(np.quantile(np.asarray(values, dtype=np.float64), quantile))


def _new_group(row: Mapping[str, str], sequence: str) -> dict[str, Any]:
    return {
        "sequence": sequence,
        "target_names": set(),
        "uniprot_ids": set(),
        "pdb_complex_ids": set(),
        "pdb_target_ids": set(),
        "dois": set(),
        "pmids": set(),
        "authors": set(),
        "publication_years": set(),
        "data_sources": set(),
        "ligands": defaultdict(list),
    }


def curate_bindingdb_series(
    *,
    bindingdb_tsv: Path,
    minimum_ligands: int = 8,
    maximum_replicate_range_pkd: float = 0.50,
    local_structure_dir: Path | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    """Stream BindingDB and return deep-series summaries, observations, and audit.

    Grouping fields are frozen in the Gate-3 protocol. Exact ligand replicates are
    collapsed to a median pKd only when their range is within the frozen tolerance.
    Every contributing BindingDB row identifier remains attached to the consensus.
    """

    if minimum_ligands < 2:
        raise CurationError("minimum_ligands must be at least two")
    counters: Counter[str] = Counter()
    groups: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    with bindingdb_tsv.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames or ())
        if missing:
            raise CurationError(f"BindingDB archive is missing columns: {sorted(missing)}")
        for row in reader:
            counters["source_rows"] += 1
            kd_nm = exact_positive_number(row["Kd (nM)"])
            if kd_nm is None:
                counters["excluded_nonexact_or_nonpositive_Kd"] += 1
                continue
            counters["exact_positive_Kd_rows"] += 1
            if (
                row[
                    "Number of Protein Chains in Target (>1 implies a multichain complex)"
                ].strip()
                != "1"
            ):
                counters["excluded_not_single_chain"] += 1
                continue
            sequence = "".join(row["BindingDB Target Chain Sequence 1"].split()).upper()
            if not sequence:
                counters["excluded_missing_target_sequence"] += 1
                continue
            publication = publication_id(row["Article DOI"], row["PMID"])
            if not publication:
                counters["excluded_missing_stable_publication"] += 1
                continue
            smiles = row["Ligand SMILES"].strip()
            if not smiles:
                counters["excluded_missing_smiles"] += 1
                continue
            try:
                canonical_smiles, scaffold = _canonical_ligand(smiles)
            except ValueError:
                counters["excluded_invalid_smiles"] += 1
                continue
            sequence_sha = hashlib.sha256(sequence.encode("ascii", "strict")).hexdigest()
            ph = normalize_condition(row["pH"])
            temperature = normalize_condition(row["Temp (C)"])
            key = (sequence_sha, publication, ph, temperature)
            group = groups.setdefault(key, _new_group(row, sequence))
            group["target_names"].add(row["Target Name"].strip())
            group["uniprot_ids"].update(
                split_accessions(
                    row["UniProt (SwissProt) Primary ID of Target Chain 1"],
                    row["UniProt (TrEMBL) Primary ID of Target Chain 1"],
                )
            )
            group["pdb_complex_ids"].update(
                split_pdb_ids(row["PDB ID(s) for Ligand-Target Complex"])
            )
            group["pdb_target_ids"].update(
                split_pdb_ids(row["PDB ID(s) of Target Chain 1"])
            )
            doi = normalize_doi(row["Article DOI"])
            if doi:
                group["dois"].add(doi)
            pmid = row["PMID"].strip()
            if pmid:
                group["pmids"].add(pmid)
            for field, name in (
                ("Authors", "authors"),
                ("Date of publication", "publication_years"),
                ("Curation/DataSource", "data_sources"),
            ):
                value = row[field].strip()
                if value:
                    group[name].add(value)
            ligand_id = row["Ligand InChI Key"].strip() or canonical_smiles
            group["ligands"][ligand_id].append(
                {
                    "reactant_set_id": row["BindingDB Reactant_set_id"].strip(),
                    "canonical_smiles": canonical_smiles,
                    "scaffold": scaffold,
                    "kd_nm": kd_nm,
                    "pkd": pkd_from_nm(kd_nm),
                }
            )
            counters["base_eligible_rows"] += 1

    local_ids: set[str] = set()
    if local_structure_dir is not None and local_structure_dir.exists():
        local_ids = {
            path.name.split(".", 1)[0].upper()
            for path in local_structure_dir.iterdir()
            if path.is_file()
        }

    summaries: list[dict[str, str]] = []
    observations: list[dict[str, str]] = []
    all_group_depths = []
    conflict_units = 0
    for key in sorted(groups):
        sequence_sha, publication, ph, temperature = key
        group = groups[key]
        retained_units = []
        for ligand_id in sorted(group["ligands"]):
            records = group["ligands"][ligand_id]
            values = [float(record["pkd"]) for record in records]
            if max(values) - min(values) > maximum_replicate_range_pkd:
                conflict_units += 1
                continue
            canonical_values = {record["canonical_smiles"] for record in records}
            if len(canonical_values) != 1:
                conflict_units += 1
                continue
            scaffolds = {record["scaffold"] for record in records}
            if len(scaffolds) != 1:
                raise CurationError("Canonical ligand has inconsistent Murcko scaffold")
            retained_units.append(
                {
                    "ligand_id": ligand_id,
                    "canonical_smiles": next(iter(canonical_values)),
                    "murcko_scaffold": next(iter(scaffolds)),
                    "pKd": median(values),
                    "replicate_count": len(records),
                    "replicate_range_pKd": max(values) - min(values),
                    "source_reactant_set_ids": ";".join(
                        sorted({record["reactant_set_id"] for record in records})
                    ),
                    "source_Kd_nM": ";".join(
                        f"{float(record['kd_nm']):.12g}" for record in records
                    ),
                }
            )
        all_group_depths.append(len(retained_units))
        if len(retained_units) < minimum_ligands:
            continue
        series_id = series_identifier(sequence_sha, publication, ph, temperature)
        pkd_values = [float(unit["pKd"]) for unit in retained_units]
        nonempty_scaffolds = {
            str(unit["murcko_scaffold"])
            for unit in retained_units
            if unit["murcko_scaffold"]
        }
        pdb_complex_ids = sorted(group["pdb_complex_ids"])
        pdb_target_ids = sorted(group["pdb_target_ids"])
        summary = {
            "series_id": series_id,
            "target_sequence_sha256": sequence_sha,
            "target_sequence_length": str(len(group["sequence"])),
            "publication_id": publication,
            "article_dois": ";".join(sorted(group["dois"])),
            "pmids": ";".join(sorted(group["pmids"])),
            "publication_years": ";".join(sorted(group["publication_years"])),
            "authors": " | ".join(sorted(group["authors"])),
            "target_names": " | ".join(sorted(group["target_names"])),
            "uniprot_ids": ";".join(sorted(group["uniprot_ids"])),
            "normalized_pH": ph,
            "normalized_temperature_C": temperature,
            "data_sources": ";".join(sorted(group["data_sources"])),
            "ligand_count": str(len(retained_units)),
            "nonempty_murcko_scaffold_count": str(len(nonempty_scaffolds)),
            "pKd_min": f"{min(pkd_values):.9g}",
            "pKd_max": f"{max(pkd_values):.9g}",
            "pKd_range": f"{max(pkd_values) - min(pkd_values):.9g}",
            "pKd_IQR": f"{_percentile(pkd_values, 0.75) - _percentile(pkd_values, 0.25):.9g}",
            "pdb_complex_ids": ";".join(pdb_complex_ids),
            "pdb_target_ids": ";".join(pdb_target_ids),
            "local_pdb_complex_ids": ";".join(sorted(set(pdb_complex_ids) & local_ids)),
            "has_holo_reference_candidate": str(bool(pdb_complex_ids)).lower(),
        }
        summaries.append(summary)
        for unit_index, unit in enumerate(
            sorted(retained_units, key=lambda item: (item["ligand_id"], item["pKd"]))
        ):
            observations.append(
                {
                    "observation_id": f"{series_id}-L{unit_index + 1:04d}",
                    "series_id": series_id,
                    "target_sequence_sha256": sequence_sha,
                    "target_sequence": group["sequence"],
                    "publication_id": publication,
                    "article_dois": summary["article_dois"],
                    "pmids": summary["pmids"],
                    "target_names": summary["target_names"],
                    "uniprot_ids": summary["uniprot_ids"],
                    "normalized_pH": ph,
                    "normalized_temperature_C": temperature,
                    "ligand_id": unit["ligand_id"],
                    "canonical_smiles": unit["canonical_smiles"],
                    "murcko_scaffold": unit["murcko_scaffold"],
                    "pKd": f"{float(unit['pKd']):.12g}",
                    "replicate_count": str(unit["replicate_count"]),
                    "replicate_range_pKd": f"{float(unit['replicate_range_pKd']):.9g}",
                    "source_reactant_set_ids": unit["source_reactant_set_ids"],
                    "source_Kd_nM": unit["source_Kd_nM"],
                    "pdb_complex_ids": summary["pdb_complex_ids"],
                    "pdb_target_ids": summary["pdb_target_ids"],
                }
            )

    summaries.sort(
        key=lambda row: (
            -int(row["has_holo_reference_candidate"] == "true"),
            -int(row["ligand_count"]),
            row["series_id"],
        )
    )
    observations.sort(key=lambda row: row["observation_id"])
    depth_counts = {
        f"series_at_least_{threshold}_ligands": sum(
            depth >= threshold for depth in all_group_depths
        )
        for threshold in (3, 5, 8, 10, 12, 15, 20, 30)
    }
    audit = {
        "schema_version": 1,
        "source": {
            "path": str(bindingdb_tsv),
            "sha256": sha256_file(bindingdb_tsv),
            "bytes": bindingdb_tsv.stat().st_size,
        },
        "contract": {
            "minimum_ligands": minimum_ligands,
            "maximum_replicate_range_pKd": maximum_replicate_range_pkd,
            "series_key": [
                "target_sequence_sha256",
                "publication_id",
                "normalized_pH",
                "normalized_temperature_C",
            ],
            "publication_is_required": True,
        },
        "row_census": dict(sorted(counters.items())),
        "series_census": {
            "all_series": len(groups),
            **depth_counts,
            "retained_deep_series": len(summaries),
            "retained_observations": len(observations),
            "retained_deep_series_with_holo_candidate": sum(
                row["has_holo_reference_candidate"] == "true" for row in summaries
            ),
            "retained_deep_series_with_local_holo_candidate": sum(
                bool(row["local_pdb_complex_ids"]) for row in summaries
            ),
            "quarantined_conflicting_ligand_units": conflict_units,
        },
    }
    return summaries, observations, audit


def select_rows(
    observations: Iterable[Mapping[str, str]], series_ids: set[str]
) -> list[dict[str, str]]:
    return [dict(row) for row in observations if row["series_id"] in series_ids]
