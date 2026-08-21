#!/usr/bin/env python3
"""Audit admitted Davis density, leakage components, ESS, and power sensitivity."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import statistics
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
import xlrd

from iscore3.gate4a.admission import connected_components
from iscore3.gate4a.labels import pkd_from_nm
from iscore3.provenance import verify_source_manifest


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _quantiles(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"minimum": None, "median": None, "maximum": None}
    return {
        "minimum": min(values),
        "median": statistics.median(values),
        "maximum": max(values),
    }


def _ligand_components(
    rows: list[dict[str, str]],
) -> tuple[list[list[str]], dict[str, str], list[tuple[str, str, str]]]:
    molecules: dict[str, Any] = {}
    scaffolds: dict[str, str] = {}
    fingerprints: dict[str, Any] = {}
    for row in rows:
        ligand_id = row["model_parent_inchikey"]
        molecule = Chem.MolFromSmiles(row["model_parent_smiles"])
        if molecule is None:
            raise RuntimeError(f"invalid admitted ligand: {row['source_name']}")
        molecules[ligand_id] = molecule
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule, includeChirality=True)
        scaffolds[ligand_id] = scaffold or f"ACYCLIC:{ligand_id[:14]}"
        fingerprints[ligand_id] = AllChem.GetMorganGenerator(
            radius=2, fpSize=2048
        ).GetFingerprint(molecule)

    ids = sorted(molecules)
    edges: list[tuple[str, str, str]] = []
    for left_index, left in enumerate(ids):
        for right in ids[left_index + 1 :]:
            reasons: list[str] = []
            if scaffolds[left] == scaffolds[right]:
                reasons.append("exact_bemis_murcko_scaffold")
            similarity = DataStructs.TanimotoSimilarity(
                fingerprints[left], fingerprints[right]
            )
            if similarity >= 0.60:
                reasons.append("morgan2_tanimoto_ge_0.60")
            if reasons:
                edges.append((left, right, ";".join(reasons)))
    components = connected_components(ids, ((left, right) for left, right, _ in edges))
    return components, scaffolds, edges


def _target_components(
    rows: list[dict[str, str]],
) -> tuple[list[list[str]], list[tuple[str, str, str]]]:
    by_id = {row["estimand_id"]: row for row in rows}
    ids = sorted(by_id)
    edges: list[tuple[str, str, str]] = []
    for left_index, left in enumerate(ids):
        for right in ids[left_index + 1 :]:
            left_row, right_row = by_id[left], by_id[right]
            reasons: list[str] = []
            family = left_row["klifs_family"]
            if family and family == right_row["klifs_family"]:
                reasons.append("same_klifs_family")
            left_pocket = left_row["klifs_pocket_sequence"]
            right_pocket = right_row["klifs_pocket_sequence"]
            identity = sum(a == b for a, b in zip(left_pocket, right_pocket)) / 85
            if identity >= 0.70:
                reasons.append("aligned_klifs85_identity_ge_0.70")
            if reasons:
                edges.append((left, right, ";".join(reasons)))
    components = connected_components(ids, ((left, right) for left, right, _ in edges))
    return components, edges


def _component_map(components: list[list[str]], prefix: str) -> dict[str, str]:
    return {
        member: f"{prefix}{index:03d}"
        for index, component in enumerate(components, start=1)
        for member in component
    }


def _write_components(
    path: Path,
    ligand_rows: list[dict[str, str]],
    receptor_rows: list[dict[str, str]],
    ligand_components: list[list[str]],
    target_components: list[list[str]],
    scaffolds: dict[str, str],
) -> None:
    ligand_map = _component_map(ligand_components, "LSC")
    target_map = _component_map(target_components, "PTC")
    ligand_sizes = Counter(ligand_map.values())
    target_sizes = Counter(target_map.values())
    output_rows: list[dict[str, Any]] = []
    for row in ligand_rows:
        entity_id = row["model_parent_inchikey"]
        component = ligand_map[entity_id]
        output_rows.append(
            {
                "entity_type": "ligand",
                "entity_id": entity_id,
                "display_name": row["source_name"],
                "component_id": component,
                "component_size": ligand_sizes[component],
                "family_or_scaffold": scaffolds[entity_id],
                "edge_policy": "Murcko equality OR Morgan2 Tanimoto >=0.60; transitive closure",
            }
        )
    for row in receptor_rows:
        entity_id = row["estimand_id"]
        component = target_map[entity_id]
        output_rows.append(
            {
                "entity_type": "receptor",
                "entity_id": entity_id,
                "display_name": row["assay_target_label"],
                "component_id": component,
                "component_size": target_sizes[component],
                "family_or_scaffold": row["klifs_family"],
                "edge_policy": (
                    "KLIFS family equality OR aligned KLIFS85 identity >=0.70; transitive closure; "
                    "structure edge pending"
                ),
            }
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite component ledger: {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compounds",
        type=Path,
        default=Path("data/processed/gate4a/davis-compound-adjudication-v1.tsv"),
    )
    parser.add_argument(
        "--receptors",
        type=Path,
        default=Path("data/processed/gate4a/davis-receptor-admission-v1.tsv"),
    )
    parser.add_argument(
        "--affinity-table",
        type=Path,
        default=Path("data/raw/gate4a/davis2011/supplementary_table_4.xls"),
    )
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=Path("data/manifests/gate4a/source-files-v1.json"),
    )
    parser.add_argument(
        "--components-output",
        type=Path,
        default=Path("data/splits/gate4a/davis-admission-components-v1.tsv"),
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=Path("reports/gate4a/evidence/davis-admitted-design-audit-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verify_source_manifest(root / args.source_manifest, repository_root=root)
    ligand_rows = [
        row for row in _read_tsv(root / args.compounds) if row["decision"] == "ACCEPTED_PARENT"
    ]
    receptor_rows = [
        row
        for row in _read_tsv(root / args.receptors)
        if row["primary_decision"] == "ACCEPTED_REFERENCE_DOMAIN"
    ]
    ligand_components, scaffolds, ligand_edges = _ligand_components(ligand_rows)
    target_components, target_edges = _target_components(receptor_rows)
    _write_components(
        root / args.components_output,
        ligand_rows,
        receptor_rows,
        ligand_components,
        target_components,
        scaffolds,
    )

    workbook = xlrd.open_workbook(str(root / args.affinity_table), on_demand=True)
    sheet = workbook.sheet_by_index(0)
    receptor_by_row = {int(row["matrix_row_1_based"]): row for row in receptor_rows}
    ligand_by_column = {
        row["affinity_matrix_name"]: row for row in ligand_rows
    }
    ligand_columns = {
        column: ligand_by_column[str(sheet.cell_value(0, column))]
        for column in range(3, sheet.ncols)
        if str(sheet.cell_value(0, column)) in ligand_by_column
    }
    exact_count = 0
    censored_count = 0
    exact_pkd: list[float] = []
    exact_by_target: Counter[str] = Counter()
    exact_by_ligand: Counter[str] = Counter()
    exact_scaffolds_by_target: dict[str, set[str]] = defaultdict(set)
    target_map = _component_map(target_components, "PTC")
    ligand_map = _component_map(ligand_components, "LSC")
    exact_cross_components: set[tuple[str, str]] = set()
    for matrix_row, receptor in receptor_by_row.items():
        for column, ligand in ligand_columns.items():
            raw = sheet.cell_value(matrix_row, column)
            if raw == "":
                censored_count += 1
                continue
            value = float(raw)
            exact_count += 1
            exact_pkd.append(pkd_from_nm(value))
            exact_by_target[receptor["estimand_id"]] += 1
            exact_by_ligand[ligand["model_parent_inchikey"]] += 1
            exact_scaffolds_by_target[receptor["estimand_id"]].add(
                ligand_map[ligand["model_parent_inchikey"]]
            )
            exact_cross_components.add(
                (
                    ligand_map[ligand["model_parent_inchikey"]],
                    target_map[receptor["estimand_id"]],
                )
            )

    ligand_component_sizes = [len(component) for component in ligand_components]
    target_component_sizes = [len(component) for component in target_components]
    effective_axis_count = min(len(ligand_components), len(target_components))
    z_sum_two_sided_alpha_0_05_power_0_80 = 1.959963984540054 + 0.8416212335729143
    mde = {
        f"paired_component_contrast_sd_{sd:g}": (
            z_sum_two_sided_alpha_0_05_power_0_80 * sd / math.sqrt(effective_axis_count)
        )
        for sd in (0.05, 0.10, 0.20, 0.30, 0.50)
    }
    exact_per_target = [exact_by_target[row["estimand_id"]] for row in receptor_rows]
    scaffolds_per_target = [
        len(exact_scaffolds_by_target[row["estimand_id"]]) for row in receptor_rows
    ]
    exact_per_ligand = [
        exact_by_ligand[row["model_parent_inchikey"]] for row in ligand_rows
    ]
    audit = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "admitted_matrix": {
            "ligand_count": len(ligand_rows),
            "receptor_count": len(receptor_rows),
            "tested_pair_count": len(ligand_rows) * len(receptor_rows),
            "exact_kd_count": exact_count,
            "right_censored_kd_count": censored_count,
            "exact_fraction": exact_count / (exact_count + censored_count),
            "exact_pkd": _quantiles(exact_pkd),
            "exact_measurements_per_target": _quantiles(
                [float(value) for value in exact_per_target]
            ),
            "exact_measurements_per_ligand": _quantiles(
                [float(value) for value in exact_per_ligand]
            ),
            "targets_with_at_least_8_exact_ligands": sum(
                value >= 8 for value in exact_per_target
            ),
            "targets_with_at_least_15_exact_ligands": sum(
                value >= 15 for value in exact_per_target
            ),
            "exact_scaffold_components_per_target": _quantiles(
                [float(value) for value in scaffolds_per_target]
            ),
            "targets_with_at_least_8_exact_scaffold_components": sum(
                value >= 8 for value in scaffolds_per_target
            ),
        },
        "ligand_leakage_components": {
            "component_count": len(ligand_components),
            "edge_count": len(ligand_edges),
            "largest_component_size": max(ligand_component_sizes),
            "singleton_component_count": sum(value == 1 for value in ligand_component_sizes),
            "policy": "Murcko equality OR Morgan2/2048 Tanimoto >=0.60; transitive closure",
        },
        "target_leakage_components": {
            "component_count_before_structure_edges": len(target_components),
            "edge_count_before_structure_edges": len(target_edges),
            "largest_component_size": max(target_component_sizes),
            "singleton_component_count": sum(value == 1 for value in target_component_sizes),
            "policy": (
                "KLIFS family equality OR aligned KLIFS85 identity >=0.70; transitive closure"
            ),
            "structure_edge_status": (
                "BLOCKED until coordinate-qualified predicted/apo pockets exist; add an edge "
                "for US-align pocket TM-score >=0.75 or Foldseek E-value <=1e-5, then recompute "
                "transitive closure without splitting any current component"
            ),
        },
        "effective_sample_size": {
            "independent_ligand_axis_upper_bound": len(ligand_components),
            "independent_target_axis_upper_bound_before_structure_edges": len(target_components),
            "conservative_two_way_axis_count": effective_axis_count,
            "observed_exact_ligand_target_component_pairs": len(exact_cross_components),
            "warning": "Cells are not independent; 23,322 tested pairs is not the ESS.",
        },
        "power_sensitivity_no_model_fitting": {
            "method": (
                "Normal-approximation two-sided paired contrast at alpha=0.05 and 80% power; "
                "n is the smaller leakage-component axis. Values are design sensitivity, not "
                "a progression threshold and not a substitute for the preregistered "
                "multiway bootstrap."
            ),
            "effective_axis_count": effective_axis_count,
            "minimum_detectable_mean_contrast_by_assumed_component_sd": mde,
        },
        "admission_decision": (
            "BLOCKED_FOR_INTERACTION until structural edges and a qualified external ledger "
            "are frozen; descriptive Davis density passes ligand-3D-only feasibility"
        ),
        "information_boundary": {
            "predictive_model_fit": False,
            "test_set_accessed": False,
            "affinity_used_for_identity_receptor_or_pocket_selection": False,
        },
    }
    output = root / args.audit_output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite admitted design audit: {output}")
    output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"Admitted Davis design: {len(ligand_rows)} ligands x {len(receptor_rows)} receptors; "
        f"{len(ligand_components)} ligand and {len(target_components)} target components"
    )


if __name__ == "__main__":
    main()
