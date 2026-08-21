#!/usr/bin/env python3
"""Freeze outcome-withheld OKL and KIRHub confirmation eligibility ledgers."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
from io import BytesIO
import json
from pathlib import Path
import re
from zipfile import ZipFile

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

from audit_kirhub_confirmation import _records_after_header
from audit_okl_confirmation import _shared_strings, _table_records


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _normalize(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def _logical_path(path: Path, root: Path) -> str:
    """Use repository-relative paths when possible and absolute replay paths otherwise."""

    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--okl-archive", type=Path, default=Path("data/raw/gate4a/confirmation/okl-supplementary-files-v2.zip"))
    parser.add_argument("--kirhub-workbook", type=Path, default=Path("data/raw/gate4a/confirmation/kirhub-supplementary-tables-v1.xlsx"))
    parser.add_argument("--davis-compounds", type=Path, default=Path("data/processed/gate4a/davis-compound-adjudication-v1.tsv"))
    parser.add_argument("--davis-receptors", type=Path, default=Path("data/processed/gate4a/davis-receptor-admission-v1.tsv"))
    parser.add_argument("--components", type=Path, default=Path("data/splits/gate4a/davis-receptor-components-final-v1.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/gate4a/confirmation"))
    parser.add_argument("--audit", type=Path, default=Path("reports/gate4a/evidence/confirmation-ledger-freeze-v1.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    davis = [row for row in _read(root / args.davis_compounds) if row["decision"] == "ACCEPTED_PARENT"]
    davis_receptors = [row for row in _read(root / args.davis_receptors) if row["primary_decision"] == "ACCEPTED_REFERENCE_DOMAIN"]
    components = {row["estimand_id"]: row["component_id"] for row in _read(root / args.components)}
    receptor_by_gene = {row["gene_symbol"]: row for row in davis_receptors}
    davis_molecules = [(row["model_parent_inchikey"], Chem.MolFromSmiles(row["model_parent_smiles"])) for row in davis]
    davis_fingerprints = {key: AllChem.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(molecule) for key, molecule in davis_molecules}
    davis_scaffolds = {key: MurckoScaffold.MurckoScaffoldSmiles(mol=molecule, includeChirality=True) or f"ACYCLIC:{key[:14]}" for key, molecule in davis_molecules}

    with ZipFile(root / args.okl_archive) as package:
        workbook_bytes = package.read("media-1.xlsx")
    with ZipFile(BytesIO(workbook_bytes)) as workbook:
        strings = _shared_strings(workbook)
        okl_compounds_raw = list(_table_records(workbook, 1, strings))
        okl_targets_raw = list(_table_records(workbook, 2, strings))
    okl_compounds = []
    for record in okl_compounds_raw:
        smiles = str(record.get("Smiles") or "")
        reported_key = str(record.get("Inchi Key") or "")
        molecule = Chem.MolFromSmiles(smiles)
        computed_key = Chem.MolToInchiKey(molecule) if molecule else ""
        exact = computed_key == reported_key and bool(reported_key)
        single = molecule is not None and len(Chem.GetMolFrags(molecule)) == 1
        fingerprint = AllChem.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(molecule) if molecule else None
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule, includeChirality=True) if molecule else ""
        exact_overlap = reported_key in davis_fingerprints
        scaffold_edge = any(scaffold and scaffold == value for value in davis_scaffolds.values())
        maximum_similarity = max((DataStructs.TanimotoSimilarity(fingerprint, value) for value in davis_fingerprints.values()), default=0.0) if fingerprint else 0.0
        leakage = exact_overlap or scaffold_edge or maximum_similarity >= 0.60
        status = "ELIGIBLE_IDENTITY_PENDING_MANUAL_QA" if exact and single and not leakage else ("EXCLUDED_DAVIS_LIGAND_LEAKAGE" if leakage else "BLOCKED_IDENTITY")
        okl_compounds.append({"panel": "OKL2026", "compound_id": str(record.get("HMSLID") or ""), "name": str(record.get("Name") or ""), "reported_smiles": smiles, "reported_inchikey": reported_key, "computed_inchikey": computed_key, "structure_roundtrip_exact": str(exact).lower(), "single_covalent_component": str(single).lower(), "exact_davis_parent_overlap": str(exact_overlap).lower(), "exact_davis_murcko_overlap": str(scaffold_edge).lower(), "maximum_davis_morgan2_tanimoto": f"{maximum_similarity:.8f}", "eligibility_status": status})
    okl_targets = []
    for record in okl_targets_raw:
        label = str(record.get("DiscoveRx Gene Symbol") or "")
        gene = str(record.get("Entrez Gene Symbol") or "")
        state = str(record.get("Wildtype or mutant") or "")
        state_specific = any(token in label.lower() for token in ("phosphorylated", "cyclin"))
        davis_row = receptor_by_gene.get(gene)
        if state != "WT" or state_specific:
            status = "EXCLUDED_NONSTANDARD_STATE_OR_MUTANT"
        elif davis_row is None:
            status = "BLOCKED_OKL_ONLY_RECEPTOR_PROVENANCE"
        else:
            status = "MAPPED_BUT_NOT_DOUBLE_COLD_FROM_DAVIS"
        okl_targets.append({"panel": "OKL2026", "assay_target_label": label, "gene_symbol": gene, "reported_state": state, "standardized_estimand_id": davis_row["estimand_id"] if davis_row else "", "davis_structure_component": components.get(davis_row["estimand_id"], "") if davis_row else "", "eligibility_status": status})

    with ZipFile(root / args.kirhub_workbook) as workbook:
        strings = _shared_strings(workbook)
        _, kir_compounds_raw = _records_after_header(workbook, strings, 1, 11)
        _, kir_targets_raw = _records_after_header(workbook, strings, 2, 21, skip_rows={22})
    davis_names = {_normalize(row["source_name"]) for row in davis} | {_normalize(row["affinity_matrix_name"]) for row in davis}
    kir_compounds = [{"panel": "KIRHub2026", "name": str(record.get("Drugs") or ""), "machine_readable_structure_present": "false", "simple_name_overlap_with_davis": str(_normalize(str(record.get("Drugs") or "")) in davis_names).lower(), "eligibility_status": "BLOCKED_NO_MACHINE_READABLE_PARENT_STRUCTURE"} for record in kir_compounds_raw]
    kir_targets = []
    accession_to_rows: dict[str, list[dict[str, str]]] = {}
    for row in davis_receptors:
        accession_to_rows.setdefault(row["uniprot_accession"], []).append(row)
    for record in kir_targets_raw:
        accession = str(record.get("Protein Accession #") or "").strip()
        matches = accession_to_rows.get(accession, [])
        construct = str(record.get("Clone") or "")
        has_range = bool(re.search(r"\baa\s*\d+\s*-\s*\d+\b", construct, re.IGNORECASE)) or "full-length" in construct.lower()
        status = "MAPPED_CONSTRUCT_PENDING_SEQUENCE_VALIDATION" if len(matches) == 1 and has_range else ("BLOCKED_NO_UNIQUE_DAVIS_ESTIMAND" if len(matches) != 1 else "BLOCKED_CONSTRUCT_BOUNDARY")
        kir_targets.append({"panel": "KIRHub2026", "assay_target_label": str(record.get("RBC Name") or ""), "gene_symbol": str(record.get("HUGO symbol") or ""), "uniprot_accession": accession, "reported_clone": construct, "standardized_estimand_id": matches[0]["estimand_id"] if len(matches) == 1 else "", "davis_structure_component": components.get(matches[0]["estimand_id"], "") if len(matches) == 1 else "", "eligibility_status": status})

    output_dir = root / args.output_dir
    paths = {
        "okl_compounds": output_dir / "okl-compound-eligibility-v1.tsv",
        "okl_targets": output_dir / "okl-target-eligibility-v1.tsv",
        "kirhub_compounds": output_dir / "kirhub-compound-eligibility-v1.tsv",
        "kirhub_targets": output_dir / "kirhub-target-eligibility-v1.tsv",
    }
    for key, records in (("okl_compounds", okl_compounds), ("okl_targets", okl_targets), ("kirhub_compounds", kir_compounds), ("kirhub_targets", kir_targets)):
        _write(paths[key], records)
    audit = {"schema_version": 1, "phase": "gate4a_provenance_closure", "decision": "LEDGERS_FROZEN_CONFIRMATION_LABELS_WITHHELD", "okl": {"compound_count": len(okl_compounds), "compound_status_counts": dict(Counter(row["eligibility_status"] for row in okl_compounds)), "target_count": len(okl_targets), "target_status_counts": dict(Counter(row["eligibility_status"] for row in okl_targets)), "strict_double_cold_pair_count": 0, "strict_double_cold_reason": "All mapped Davis targets fall within the dominant structure-leakage component; OKL-only targets lack closed receptor provenance. No pair is released."}, "kirhub": {"compound_count": len(kir_compounds), "compound_status_counts": dict(Counter(row["eligibility_status"] for row in kir_compounds)), "target_count": len(kir_targets), "target_status_counts": dict(Counter(row["eligibility_status"] for row in kir_targets)), "eligible_pair_count": 0, "reason": "Compound structures are absent from the source supplement."}, "information_boundary": {"outcome_tables_read": False, "cell_level_outcomes_written": False, "used_for_tuning": False, "labels_released": False}, "ledgers": {key: {"path": _logical_path(path, root), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "rows": sum(1 for _ in path.open()) - 1} for key, path in paths.items()}}
    audit_path = root / args.audit
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    if audit_path.exists():
        raise FileExistsError(f"refusing to overwrite {audit_path}")
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
